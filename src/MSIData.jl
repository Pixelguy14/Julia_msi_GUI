# src/MSIData.jl

"""
This file defines the unified MSIData object and the associated data access layer,
including caching and iteration logic, for handling large mzML and imzML datasets
efficiently.
"""

using Base64, Libz, Serialization, Printf, DataFrames, Base.Threads, StatsBase, Mmap

const FILE_HANDLE_LOCK = ReentrantLock()


# Abstract type for different data sources (e.g., mzML, imzML)
# This allows dispatching to the correct binary reading logic.
abstract type MSDataSource end

# Concrete type for imzML data sources
"""
    ImzMLSource <: MSDataSource

A data source for `.imzML` files, holding a handle to the binary `.ibd` file
and the expected format for m/z and intensity arrays.
"""
struct ImzMLSource <: MSDataSource
    ibd_handles::Vector{IO} # HandlePool: One handle per thread
    mz_format::Type
    intensity_format::Type
    mmap_data::Union{Vector{UInt8}, Nothing}
    is_any_compressed::Bool # Cached for zero-allocation dispatch
end

"""
    MzMLSource <: MSDataSource

A data source for `.mzML` files, holding a handle to the `.mzML` file itself
(which contains the binary data encoded in Base64) and the expected data formats.
"""
struct MzMLSource <: MSDataSource
    file_handles::Vector{IO} # HandlePool: One handle per thread
    mz_format::Type
    intensity_format::Type
    mmap_data::Union{Vector{UInt8}, Nothing}
end

# --- HandlePool Helpers --- #
"""
    get_handle(source::ImzMLSource) -> IO
    get_handle(source::MzMLSource) -> IO

Retrieves a thread-local file handle from the source's pool.
"""
function get_handle(source::ImzMLSource)
    tid = Threads.threadid()
    if tid <= length(source.ibd_handles)
        return source.ibd_handles[tid]
    else
        return source.ibd_handles[1]
    end
end

function get_handle(source::MzMLSource)
    tid = Threads.threadid()
    if tid <= length(source.file_handles)
        return source.file_handles[tid]
    else
        return source.file_handles[1]
    end
end

"""
    SpectrumAsset

Contains metadata for a single binary data array (m/z or intensity) within a spectrum.

# Fields
- `format`: The data type of the elements (e.g., `Float32`, `Int64`).
- `is_compressed`: A boolean flag indicating if the data is compressed (e.g., with zlib).
- `offset`: The byte offset of the data within the file (`.ibd` for imzML, `.mzML` for mzML).
- `encoded_length`: The length of the data, which has different meanings depending on the format:
    - **Uncompressed imzML**: The number of elements (points) in the array.
    - **Compressed imzML**: The number of bytes of the compressed data chunk in the `.ibd` file.
    - **mzML**: The number of characters in the Base64 encoded string within the `.mzML` file.
- `axis_type`: A symbol (`:mz` or `:intensity`) indicating the type of data.
"""
struct SpectrumAsset
    format::Type
    is_compressed::Bool
    offset::Int64
    encoded_length::Int32
    # For mzML, axis_type is needed to distinguish mz from intensity.
    # For imzML, this can be ignored as the order is fixed.
    axis_type::Symbol
    
    # Pre-computed analytics
    min_val::Float64
    max_val::Float64
end

"""
    SpectrumMode

An enum representing the type of a spectrum's data.

- `CENTROID`: The spectrum contains centroided data (i.e., processed peaks).
- `PROFILE`: The spectrum contains raw profile data.
- `UNKNOWN`: The spectrum type is not specified.
"""
@enum SpectrumMode CENTROID=1 PROFILE=2 UNKNOWN=3

"""
    InstrumentMetadata

Contains metadata about the instrument configuration and acquisition parameters.
This information is parsed from the file header and is crucial for guiding
preprocessing steps.

# Fields
- `resolution`: Instrument resolution at a specific m/z, if available.
- `acquisition_mode`: The overall data acquisition mode (`:profile`, `:centroid`, or `:mixed`).
- `mz_axis_type`: The nature of the m/z axis (`:regular`, `:irregular`, `:mixed`).
- `calibration_status`: The calibration state of the data (`:internal`, `:external`, `:uncalibrated`).
- `instrument_model`: The vendor and model of the instrument.
- `mass_accuracy_ppm`: Manufacturer-specified mass accuracy in parts-per-million.
- `laser_settings`: A dictionary of MSI-specific laser parameters.
- `polarity`: The polarity of the acquisition (`:positive`, `:negative`, `:unknown`).
"""
struct InstrumentMetadata
    resolution::Union{Float64, Nothing}
    acquisition_mode::Symbol
    mz_axis_type::Symbol
    calibration_status::Symbol
    instrument_model::String
    mass_accuracy_ppm::Union{Float64, Nothing}
    laser_settings::Union{Dict, Nothing}
    polarity::Symbol
    vendor_preprocessing_steps::Union{Vector{String}, Nothing} # New field for provenance
end

# Default constructor for initializing with unknown values
function InstrumentMetadata()
    return InstrumentMetadata(
        nothing,        # resolution
        :unknown,       # acquisition_mode
        :unknown,       # mz_axis_type
        :uncalibrated,  # calibration_status
        "Not Available",# instrument_model
        nothing,        # mass_accuracy_ppm
        nothing,        # laser_settings
        :unknown,       # polarity
        nothing         # vendor_preprocessing
    )
end

"""
    SpectrumMetadata

Contains all metadata for a single spectrum, common to both imzML and mzML formats.

# Fields
- `x`, `y`: The spatial coordinates of the spectrum (for imzML only).
- `id`: The unique identifier string for the spectrum (for mzML only).
- `mode`: The spectrum mode (`CENTROID` or `PROFILE`).
- `mz_asset`: A `SpectrumAsset` for the m/z array.
- `int_asset`: A `SpectrumAsset` for the intensity array.
"""
struct SpectrumMetadata
    # For imzML
    x::Int32
    y::Int32
    
    # For mzML
    id::String
    type::Symbol # NEW: e.g., :sample, :blank, :qc

    mode::SpectrumMode
    # Common binary data info
    mz_asset::SpectrumAsset
    int_asset::SpectrumAsset
end

"""
    SpectrumMetadataBinary

A fixed-size version of SpectrumMetadata for high-speed binary serialization.
Used for the metadata cache (.cache files).
"""
struct SpectrumMetadataBinary
    x::Int32
    y::Int32
    mode::Int8
    mz_offset::Int64
    mz_encoded_len::Int32
    int_offset::Int64
    int_encoded_len::Int32
    min_mz::Float32 # Persistent analytics
    max_mz::Float32 # Persistent analytics
end

"""
    MSIData

The primary object for interacting with mass spectrometry data. It provides a unified
interface for both `.mzML` and `.imzML` files and includes an LRU cache for
efficient repeated access to spectra.

# Fields
- `source`: The underlying `MSDataSource` (`ImzMLSource` or `MzMLSource`).
- `spectra_metadata`: A vector of `SpectrumMetadata` for all spectra in the file.
- `image_dims`: A tuple `(width, height)` of the spatial dimensions (for imzML).
- `coordinate_map`: A matrix mapping `(x, y)` coordinates to a linear spectrum index (for imzML).
- `cache`: A dictionary holding cached spectra, mapping index to `(mz, intensity)`.
- `cache_order`: A vector of indices tracking usage for the LRU cache policy.
- `cache_size`: The maximum number of spectra to store in the cache.
- `cache_lock`: A `ReentrantLock` to ensure thread-safe access to the cache.
- `global_min_mz`: Cached global minimum m/z value across all spectra.
- `global_max_mz`: Cached global maximum m/z value across all spectra.
- `spectrum_stats_df`: A `DataFrame` containing pre-computed per-spectrum analytics (e.g., TIC, BPI).
"""
mutable struct MSIData
    source::MSDataSource
    spectra_metadata::Vector{SpectrumMetadata}
    instrument_metadata::Union{InstrumentMetadata, Nothing}
    image_dims::Tuple{Int, Int}
    coordinate_map::Union{Matrix{Int}, Nothing}

    # LRU Cache for GetSpectrum
    cache::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}
    cache_order::Vector{Int}
    cache_size::Int
    cache_lock::ReentrantLock

    # Buffer Pool for binary data operations
    buffer_pool::SimpleBufferPool # Existing UInt8 pool (mostly for mzML base64)

    # Pre-computed analytics/metadata - use Base.Threads.Atomic for compatibility
    global_min_mz::Base.Threads.Atomic{Float64}
    global_max_mz::Base.Threads.Atomic{Float64}
    spectrum_stats_df::Union{DataFrame, Nothing}
    bloom_filters::Union{Vector{Union{BloomFilter{Int}, Nothing}}, Nothing}
    analytics_ready::AtomicFlag
    preprocessing_hints::Union{Dict{Symbol, Any}, Nothing} # For auto-determined parameters

    function MSIData(source, metadata, instrument_meta, dims, coordinate_map, cache_size)
        obj = new(source, metadata, instrument_meta, dims, coordinate_map, 
                  Dict(), [], min(10, cache_size), ReentrantLock(),
                  SimpleBufferPool(),
                  Base.Threads.Atomic{Float64}(0.0), Base.Threads.Atomic{Float64}(0.0), 
                  nothing, nothing, AtomicFlag(), nothing)
        
        # Initialize mz bounds cleanly instead of Inf
        Base.Threads.atomic_xchg!(obj.global_min_mz, 1e9)
        Base.Threads.atomic_xchg!(obj.global_max_mz, -1e9)

        
        # Ensure all file handles in the pool are closed when the object is garbage collected
        finalizer(obj) do o
            if o.source isa ImzMLSource
                for h in o.source.ibd_handles
                    isopen(h) && close(h)
                end
            elseif o.source isa MzMLSource
                for h in o.source.file_handles
                    isopen(h) && close(h)
                end
            end
        end
        return obj
    end
end

"""
    set_global_mz_range!(data::MSIData, min_mz::Float64, max_mz::Float64)

Sets the global m/z range for the MSIData object.

# Arguments

- `data::MSIData`: The MSIData object.
- `min_mz::Float64`: The minimum m/z value.
- `max_mz::Float64`: The maximum m/z value.

# Returns
- `nothing`
"""
function set_global_mz_range!(data::MSIData, min_mz::Float64, max_mz::Float64)
    Base.Threads.atomic_xchg!(data.global_min_mz, min_mz)
    Base.Threads.atomic_xchg!(data.global_max_mz, max_mz)
end

"""
    get_global_mz_range(data::MSIData)

Gets the global m/z range for the MSIData object.

# Arguments

- `data::MSIData`: The MSIData object.

# Returns
- `(min_mz, max_mz)`: The global m/z range.
"""
function get_global_mz_range(data::MSIData)
    # Atomic read using atomic_add! with zero
    min_mz = Base.Threads.atomic_add!(data.global_min_mz, 0.0)
    max_mz = Base.Threads.atomic_add!(data.global_max_mz, 0.0)
    return (min_mz, max_mz)
end

"""
    set_analytics_data!(data::MSIData, stats_df::DataFrame, bloom_filters::Vector)

Sets the analytics data for the MSIData object.

# Arguments

- `data::MSIData`: The MSIData object.
- `stats_df::DataFrame`: The statistics DataFrame.
- `bloom_filters::Vector`: The bloom filters.

# Returns
- `nothing`
"""
function set_analytics_data!(data::MSIData, stats_df::DataFrame, bloom_filters::Vector)
    lock(data.cache_lock) do
        data.spectrum_stats_df = stats_df
        data.bloom_filters = bloom_filters
    end
end

"""
    get_bloom_filters(data::MSIData)

Gets the bloom filters for the MSIData object.

# Arguments

- `data::MSIData`: The MSIData object.

# Returns
- `bloom_filters::Vector`: The bloom filters.
"""
function get_bloom_filters(data::MSIData)
    lock(data.cache_lock) do
        return data.bloom_filters
    end
end

"""
    get_spectrum_stats(data::MSIData)

Gets the spectrum statistics for the MSIData object.

# Arguments

- `data::MSIData`: The MSIData object.

# Returns
- `stats_df::DataFrame`: The statistics DataFrame.
"""
function get_spectrum_stats(data::MSIData)
    return data.spectrum_stats_df
end

"""
    Base.close(data::MSIData)

Explicitly closes the file handles associated with the `MSIData` object and clears the spectrum cache.
It is good practice to call this method when you are finished with an `MSIData` object to release resources immediately.

# Arguments

- `data::MSIData`: The MSIData object.

# Returns
- `nothing`
"""
function Base.close(data::MSIData)
    if data.source isa ImzMLSource
        for handle in data.source.ibd_handles
            isopen(handle) && close(handle)
        end
    elseif data.source isa MzMLSource
        for handle in data.source.file_handles
            isopen(handle) && close(handle)
        end
    end

    # Clear cache
    empty!(data.cache)
    empty!(data.cache_order)

    println("MSIData object closed and resources released.")
end



"""
    get_masked_spectrum_indices(data::MSIData, mask_matrix::BitMatrix) -> Set{Int}

Converts a 2D boolean mask matrix (e.g., loaded from a PNG) into a `Set` of linear
spectrum indices that fall within the `true` regions of the mask.

# Arguments
- `data`: The `MSIData` object containing the `coordinate_map`.
- `mask_matrix`: A `BitMatrix` where `true` indicates a pixel is part of the ROI.

# Returns
- A `Set{Int}` containing the linear indices of spectra within the mask.
"""
function get_masked_spectrum_indices(data::MSIData, mask_matrix::BitMatrix)
    if data.coordinate_map === nothing
        throw(ArgumentError("Coordinate map not available. Cannot apply mask to non-imaging data."))
    end
    
    width, height = data.image_dims
    mask_height, mask_width = size(mask_matrix)

    if mask_width != width || mask_height != height
        throw(ArgumentError("Mask dimensions ($(mask_width)x$(mask_height)) do not match image dimensions ($(width)x$(height))."))
    end

    masked_indices = Set{Int}()
    for y in 1:height
        for x in 1:width
            if mask_matrix[y, x] # If this pixel is part of the mask
                idx = data.coordinate_map[x, y]
                if idx != 0 # Ensure there's an actual spectrum at this coordinate
                    push!(masked_indices, idx)
                end
            end
        end
    end
    return masked_indices
end


# --- Internal function for reading binary data --- #

"""
    read_binary_vector(io::IO, asset::SpectrumAsset)

Reads and decodes a single binary data vector (e.g., m/z or intensity array)
from a `.mzML` file. The data is expected to be Base64-encoded and may be
compressed.

This internal function handles:
1. Reading the raw Base64 string from the file at the specified offset.
2. Decoding the Base64 string into bytes.
3. Decompressing the bytes using zlib if `asset.is_compressed` is true.
4. Interpreting the resulting bytes as a vector of the specified format.

Note: Byte order conversion (e.g., from little-endian to host) is not performed
by this function and is assumed to be handled by the caller if necessary.

# Arguments
- `io`: The IO stream of the `.mzML` file.
- `asset`: The `SpectrumAsset` containing metadata for the array.

# Returns
- A `Vector` of the appropriate type containing the decoded data.
"""
function read_binary_vector(data::MSIData, io::IO, asset::SpectrumAsset)
    if asset.offset < 0
        throw(FileFormatError("Invalid asset offset: $(asset.offset)"))
    end
    
    # Use mmap view if available for Base64 (mzML)
    b64_string = (data.source isa MzMLSource && data.source.mmap_data !== nothing) ? 
                 String(view(data.source.mmap_data, (asset.offset + 1):(asset.offset + asset.encoded_length))) :
                 String(read(seek(io, asset.offset), asset.encoded_length))
    
    local decoded_bytes::Vector{UInt8}
    
    if asset.is_compressed
        temp_decoded = Base64.base64decode(b64_string)
        decoded_bytes = Libz.inflate(temp_decoded)
    else
        decoded_bytes = Base64.base64decode(b64_string)
    end
    
    alignment = sizeof(asset.format)
    n_elements = length(decoded_bytes) ÷ alignment

    if n_elements * alignment != length(decoded_bytes)
        throw(FileFormatError("Size of decoded byte array is not a multiple of the element size."))
    end
    
    reinterpreted_array = reinterpret(asset.format, decoded_bytes)
    
    # Optimization: Use a temporary array for byte order conversion.
    # We could use the ResourcePool here if we wanted to return a Float64 vector,
    # but currently we return the native format.
    out_array = [ltoh(x) for x in reinterpreted_array]
    
    return out_array
end



# Overload for different source types
"""
    read_spectrum_from_disk(source::ImzMLSource, meta::SpectrumMetadata)

Reads a single spectrum's m/z and intensity arrays directly from the `.ibd`
binary file for an `.imzML` dataset. It handles reading the raw binary data
and converting it from little-endian to the host's native byte order.

# Arguments
- `source`: An `ImzMLSource` containing the file handle and data formats.
- `meta`: The `SpectrumMetadata` for the spectrum to be read.

# Returns
- A tuple `(mz, intensity)` containing the two requested data arrays.
"""
@inline function read_spectrum_from_disk(source::ImzMLSource, meta::SpectrumMetadata)
    # 1. Use Mmap logic if available (implemented in previous step)
    if source.mmap_data !== nothing
        mz_offset = meta.mz_asset.offset
        mz_byte_len = sizeof(source.mz_format) * meta.mz_asset.encoded_length
        int_offset = meta.int_asset.offset
        int_byte_len = sizeof(source.intensity_format) * meta.int_asset.encoded_length

        mz_view_raw = view(source.mmap_data, (mz_offset + 1):(mz_offset + mz_byte_len))
        int_view_raw = view(source.mmap_data, (int_offset + 1):(int_offset + int_byte_len))

        mz_reinterpreted = reinterpret(source.mz_format, mz_view_raw)
        int_reinterpreted = reinterpret(source.intensity_format, int_view_raw)

        # TRUE Zero-copy logic: avoid allocations if host matches file endianness (LE for .ibd)
        if Base.ENDIAN_BOM == 0x04030201 # Little Endian Host (Common for Linux/X86)
            mz = mz_reinterpreted
            intensity = int_reinterpreted
        else
            # On Big Endian hosts, we MUST allocate and byte-swap
            mz = ltoh.(mz_reinterpreted)
            intensity = ltoh.(int_reinterpreted)
        end

        validate_spectrum_data(mz, intensity, meta.id)
        return mz, intensity
    end

    # 2. Use HandlePool logic if Mmap is not available
    handle = get_handle(source)
    
    mz = Array{source.mz_format}(undef, meta.mz_asset.encoded_length)
    intensity = Array{source.intensity_format}(undef, meta.int_asset.encoded_length)

    # Note: No lock() required here because we are using a thread-local handle!
    seek(handle, meta.mz_asset.offset)
    read!(handle, mz)
    
    seek(handle, meta.int_asset.offset)
    read!(handle, intensity)

    mz .= ltoh.(mz)
    intensity .= ltoh.(intensity)

    validate_spectrum_data(mz, intensity, meta.id) 
    return mz, intensity
end

"""
    read_spectrum_from_disk(source::MzMLSource, meta::SpectrumMetadata)

Reads a single spectrum's m/z and intensity arrays from a `.mzML` file.
This is a high-level function that orchestrates the reading process by calling
`read_binary_vector` for each of the m/z and intensity arrays.

# Arguments
- `source`: A `MzMLSource` containing the file handle.
- `meta`: The `SpectrumMetadata` for the spectrum to be read.

# Returns
- A tuple `(mz, intensity)` containing the two requested data arrays.
"""
function read_spectrum_from_disk(data::MSIData, source::MzMLSource, meta::SpectrumMetadata)
    handle = get_handle(source)
    # For mzML, data is Base64 encoded within the XML
    mz = read_binary_vector(data, handle, meta.mz_asset)
    intensity = read_binary_vector(data, handle, meta.int_asset)
    
    mz_f64, intensity_f64 = Float64.(mz), Float64.(intensity)

    validate_spectrum_data(mz, intensity, meta.id) 

    return mz, intensity
end

# --- Metadata Caching (Sprint 1: Milestone 4) --- #

"""
    save_metadata_cache(data::MSIData, cache_path::String)

Serializes the spectrum metadata to a custom binary format for near-instant loading.
"""
function save_metadata_cache(data::MSIData, cache_path::String)
    open(cache_path, "w") do io
        # Write magic number and version (v2 adds min_mz/max_mz to SpectrumMetadataBinary)
        write(io, "JMSI")
        write(io, Int32(2))
        
        # Write number of spectra
        num_spectra = length(data.spectra_metadata)
        write(io, Int32(num_spectra))
        
        # Write global formats (assuming uniform for now)
        # We'll write the names of the types as strings for safety
        write(io, string(data.source.mz_format))
        write(io, "\n")
        write(io, string(data.source.intensity_format))
        write(io, "\n")

        # Convert to binary structs and write in one block
        binary_metadata = Vector{SpectrumMetadataBinary}(undef, num_spectra)
        for i in 1:num_spectra
            m = data.spectra_metadata[i]
            binary_metadata[i] = SpectrumMetadataBinary(
                m.x, m.y, Int8(m.mode),
                m.mz_asset.offset, m.mz_asset.encoded_length,
                m.int_asset.offset, m.int_asset.encoded_length,
                Float32(m.mz_asset.min_val), Float32(m.mz_asset.max_val)
            )
        end
        write(io, binary_metadata)
    end
    @debug "Metadata cache saved to $cache_path"
end

"""
    load_metadata_cache(cache_path::String, mz_format::Type, int_format::Type) -> Vector{SpectrumMetadata}

Loads spectrum metadata from a custom binary cache file.
"""
function load_metadata_cache(cache_path::String, mz_format::Type, int_format::Type)
    open(cache_path, "r") do io
        magic = read(io, 4)
        if String(magic) != "JMSI"
            error("Invalid cache file format.")
        end
        version = read(io, Int32)
        if version != 2
            error("Unsupported cache version $version (expected 2). Delete the .cache file to regenerate.")
        end
        num_spectra = read(io, Int32)
        
        # Skip format strings (we already have them from the header or caller)
        readline(io)
        readline(io)
        
        # Read all binary metadata in one swoop
        binary_metadata = Vector{SpectrumMetadataBinary}(undef, num_spectra)
        read!(io, binary_metadata)
        
        # Convert back to SpectrumMetadata
        spectra_metadata = Vector{SpectrumMetadata}(undef, num_spectra)
        for i in 1:num_spectra
            b = binary_metadata[i]
            mz_asset = SpectrumAsset(mz_format, false, b.mz_offset, b.mz_encoded_len, :mz, Float64(b.min_mz), Float64(b.max_mz))
            int_asset = SpectrumAsset(int_format, false, b.int_offset, b.int_encoded_len, :intensity, 0.0, 0.0)
            
            spectra_metadata[i] = SpectrumMetadata(
                b.x, b.y, "", :sample, SpectrumMode(b.mode),
                mz_asset, int_asset
            )
        end
        return spectra_metadata
    end
end

# --- Public API --- #

"""
    GetSpectrum(data::MSIData, index::Int)

Retrieves a single spectrum by its index, utilizing a thread-safe LRU cache for performance.

If the spectrum is not in the cache, it is read from disk, and the cache is updated.
This function is the core of the "Indexed" and "Cache" access patterns.

# Arguments
- `data`: The `MSIData` object.
- `index`: The linear index of the spectrum to retrieve.

# Returns
- A tuple `(mz, intensity)` containing the spectrum's data arrays.
"""
function GetSpectrum(data::MSIData, index::Int)
    if index < 1 || index > length(data.spectra_metadata)
        throw(SpectrumNotFoundError(index))
    end

    # Phase 1: Check the cache (with lock)
    lock(data.cache_lock)
    try
        if haskey(data.cache, index)
            # Cache Hit: Move item to the end of the LRU list and return from cache
            filter!(x -> x != index, data.cache_order)
            push!(data.cache_order, index)
            return data.cache[index]
        end
    finally
        unlock(data.cache_lock)
    end

    # Phase 2: Cache Miss - Read from disk (no lock)
    meta = data.spectra_metadata[index]
    
    local spectrum
    if data.source isa ImzMLSource
        spectrum = read_spectrum_from_disk(data.source, meta)
    else # MzMLSource
        spectrum = read_spectrum_from_disk(data, data.source, meta)
    end

    # Phase 3: Update cache (with lock) and get final value
    return lock(data.cache_lock) do
        if haskey(data.cache, index)
            # Another thread got here first, use its result
            return data.cache[index]
        end

        # This thread is first to update cache
        if data.cache_size > 0
            if length(data.cache) >= data.cache_size
                # Evict the least recently used item (at the front of the list)
                lru_index = popfirst!(data.cache_order)
                delete!(data.cache, lru_index)
            end
            data.cache[index] = spectrum
            push!(data.cache_order, index)
        end
        return spectrum
    end
end

"""
    GetSpectrum(data::MSIData, x::Int, y::Int)

Retrieves a single spectrum by its (x, y) coordinates for imaging data (`.imzML`).
This method uses the `coordinate_map` for efficient index lookup and then calls
the indexed `GetSpectrum` method, benefiting from caching.

# Arguments
- `data`: The `MSIData` object.
- `x`: The x-coordinate of the spectrum.
- `y`: The y-coordinate of the spectrum.

# Returns
- A tuple `(mz, intensity)` containing the spectrum's data arrays.
"""
function GetSpectrum(data::MSIData, x::Int, y::Int)
    if data.coordinate_map === nothing
        throw(ArgumentError("Coordinate map not available. This method is only for imaging data loaded from .imzML files."))
    end
    width, height = data.image_dims
    if x < 1 || x > width || y < 1 || y > height
        throw(ArgumentError("Coordinates ($x, $y) out of bounds for image dimensions ($width, $height)."))
    end

    index = data.coordinate_map[x, y]
    if index == 0
        throw(SpectrumNotFoundError((x, y)))
    end

    return GetSpectrum(data, index) # Call the existing method
end

"""
    process_spectrum(f::Function, data::MSIData, index::Int)

A "function barrier" helper for safely processing a single spectrum.

This function retrieves a spectrum and then immediately calls the provided function `f`
with the resulting `(mz, intensity)` arrays. While the `GetSpectrum` API is now type-stable,
this pattern remains good practice for separating data access from data processing.

Use this function when you need to perform performance-critical operations on a
single spectrum. Your logic should be inside the function passed to this helper.

For bulk processing of all spectra, consider using the `_iterate_spectra_fast`
function, which is optimized for sequential access and avoids cache overhead.

# Arguments

- `f::Function`: The function to apply to the spectrum.
- `data::MSIData`: The MSIData object.
- `index::Int`: The index of the spectrum.

# Returns
- `result::Any`: The result of applying the function to the spectrum.

# Example

```julia
process_spectrum(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function process_spectrum(f::Function, data::MSIData, index::Int)
    # This call is type-unstable.
    mz, intensity = GetSpectrum(data, index)
    
    # By immediately calling `f`, we create a "function barrier".
    # Julia will compile a specialized version of `f` for the
    # concrete types of `mz` and `intensity`.
    return f(mz, intensity)
end

"""
    process_spectrum(f::Function, data::MSIData, x::Int, y::Int)

A coordinate-based version of the "function barrier" helper. See `process_spectrum`
for details on the performance pattern.

# Arguments

- `f::Function`: The function to apply to the spectrum.
- `data::MSIData`: The MSIData object.
- `x::Int`: The x-coordinate of the spectrum.
- `y::Int`: The y-coordinate of the spectrum.

# Returns
- `result::Any`: The result of applying the function to the spectrum.

# Example

```julia
process_spectrum(data, 1, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function process_spectrum(f::Function, data::MSIData, x::Int, y::Int)
    mz, intensity = GetSpectrum(data, x, y)
    return f(mz, intensity)
end

"""
    precompute_analytics(msi_data::MSIData)

Performs a single pass over the entire dataset to pre-compute and cache important
analytics. This function populates the `global_min_mz`, `global_max_mz`, and
`spectrum_stats_df` fields of the `MSIData` object.

The computed statistics include:
- Global minimum and maximum m/z values.
- Per-spectrum:
    - Total Ion Count (TIC)
    - Base Peak Intensity (BPI)
    - m/z of the base peak
    - Number of data points
    - Minimum and maximum m/z

Subsequent calls to functions like `get_total_spectrum` will be much faster
as they can use this cached data. This function modifies the `MSIData` object in-place
and is idempotent.

# Arguments

- `msi_data::MSIData`: The MSIData object.

# Returns
- `nothing`
"""
function precompute_analytics(msi_data::MSIData)
    # Double-checked locking pattern
    if is_set(msi_data.analytics_ready)
        return
    end
    
    lock(msi_data.cache_lock) do
        if is_set(msi_data.analytics_ready)
            return
        end
        
        println("Pre-computing analytics (streaming mode)...")
        start_time = time_ns()

        num_spectra = length(msi_data.spectra_metadata)
        
        # Initialize thread-local variables for global stats
        thread_local_min_mz = [Inf for _ in 1:Base.Threads.nthreads()]
        thread_local_max_mz = [-Inf for _ in 1:Base.Threads.nthreads()]

        # Initialize vectors for per-spectrum stats
        tics = Vector{Float64}(undef, num_spectra)
        bpis = Vector{Float64}(undef, num_spectra)
        bp_mzs = Vector{Float64}(undef, num_spectra)
        num_points = Vector{Int}(undef, num_spectra)
        min_mzs = Vector{Float64}(undef, num_spectra)
        max_mzs = Vector{Float64}(undef, num_spectra)
        modes = Vector{SpectrumMode}(undef, num_spectra)
        
        # Precompute Bloom filters during this pass
        bloom_filters = Vector{Union{BloomFilter{Int}, Nothing}}(undef, num_spectra)
        fill!(bloom_filters, nothing)

        # Process in chunks to reduce memory pressure
        chunk_size = min(10000, num_spectra)
        
        for chunk_start in 1:chunk_size:num_spectra
            chunk_end = min(chunk_start + chunk_size - 1, num_spectra)
            chunk_range = chunk_start:chunk_end
            
            println("Processing chunk $chunk_start - $chunk_end / $num_spectra")
            
            # Process current chunk
            _iterate_spectra_fast(msi_data, chunk_range) do idx, mz, intensity
                # Store metadata
                modes[idx] = msi_data.spectra_metadata[idx].mode
                
                if isempty(mz)
                    tics[idx] = 0.0
                    bpis[idx] = 0.0
                    bp_mzs[idx] = 0.0
                    num_points[idx] = 0
                    min_mzs[idx] = Inf
                    max_mzs[idx] = -Inf
                    return
                end

                # Update thread-local global m/z range
                local_min, local_max = extrema(mz)
                thread_id = Base.Threads.threadid()
                thread_local_min_mz[thread_id] = min(thread_local_min_mz[thread_id], local_min)
                thread_local_max_mz[thread_id] = max(thread_local_max_mz[thread_id], local_max)
                min_mzs[idx] = local_min
                max_mzs[idx] = local_max

                # Calculate per-spectrum stats
                tics[idx] = sum(intensity)
                max_int, max_idx = findmax(intensity)
                bpis[idx] = max_int
                bp_mzs[idx] = mz[max_idx]
                num_points[idx] = length(mz)
                bloom_filters[idx] = create_bloom_filter_for_spectrum(mz)
            end
            
            # Force GC between chunks to control memory
            if chunk_end % 5000 == 0
                GC.gc()
            end
        end

        # Combine thread-local global stats
        g_min_mz = minimum(thread_local_min_mz)
        g_max_mz = maximum(thread_local_max_mz)

        # Only update if we found valid ranges
        if isfinite(g_min_mz) && isfinite(g_max_mz) && g_min_mz < g_max_mz
            Base.Threads.atomic_xchg!(msi_data.global_min_mz, g_min_mz)
            Base.Threads.atomic_xchg!(msi_data.global_max_mz, g_max_mz)
        else
            # Fallback: calculate from first spectrum
            try
                if num_spectra > 0
                    mz, _ = GetSpectrum(msi_data, 1)
                    if !isempty(mz)
                        Base.Threads.atomic_xchg!(msi_data.global_min_mz, minimum(mz))
                        Base.Threads.atomic_xchg!(msi_data.global_max_mz, maximum(mz))
                    end
                end
            catch e
                @warn "Could not determine global m/z range from first spectrum: $e"
            end
        end
        
        # Create results
        computed_stats_df = DataFrame(
            SpectrumID = 1:num_spectra,
            TIC = tics,
            BPI = bpis,
            BasePeakMZ = bp_mzs,
            NumPoints = num_points,
            MinMZ = min_mzs,
            MaxMZ = max_mzs,
            Mode = modes
        )
        
        # Set other data
        msi_data.spectrum_stats_df = computed_stats_df
        msi_data.bloom_filters = bloom_filters  # Will be created on-demand
        
        set!(msi_data.analytics_ready)

        duration = (time_ns() - start_time) / 1e9
        @printf "Streaming analytics complete in %.2f seconds.\n" duration
    end
    
    return
end




"""
    get_total_spectrum_core(msi_data::MSIData; num_bins::Int=2000, masked_indices::Union{Set{Int}, Nothing}=nothing)

Internal function to calculate the total spectrum.

It uses a two-pass approach:
1. The first pass finds the global m/z range by iterating through all spectra.
2. The second pass sums intensities into a pre-defined number of bins.

This function is called by `get_total_spectrum`.

# Arguments

- `msi_data::MSIData`: The MSIData object.
- `num_bins::Int`: The number of bins to use for the total spectrum.
- `masked_indices::Union{Set{Int}, Nothing}`: The indices of the spectra to mask.

# Returns
- `(mz_range, total_spectrum, num_spectra_processed)`: A tuple containing the m/z range, total spectrum, and spectrum count.
"""
function get_total_spectrum_core(msi_data::MSIData; num_bins::Int=2000, masked_indices::Union{Set{Int}, Nothing}=nothing)
    println("Calculating total spectrum (2-pass method)...")
    total_start_time = time_ns()

    local global_min_mz, global_max_mz

    if Base.Threads.atomic_add!(msi_data.global_min_mz, 0.0) !== Inf
        println("  Using pre-computed m/z range.")
        global_min_mz = Base.Threads.atomic_add!(msi_data.global_min_mz, 0.0)
        global_max_mz = Base.Threads.atomic_add!(msi_data.global_max_mz, 0.0)
    else
        # 1. First Pass: Find the global m/z range
        pass1_start_time = time_ns()
        println("  Pass 1: Finding global m/z range...")
        g_min_mz = Inf
        g_max_mz = -Inf
        _iterate_spectra_fast(msi_data, masked_indices === nothing ? nothing : collect(masked_indices)) do idx, mz, _
            if !isempty(mz)
                local_min, local_max = extrema(mz)
                g_min_mz = min(g_min_mz, local_min)
                g_max_mz = max(g_max_mz, local_max)
            end
        end
        pass1_duration = (time_ns() - pass1_start_time) / 1e9
        
        if !isfinite(g_min_mz)
            @warn "Could not determine a valid m/z range. All spectra might be empty."
            return (Float64[], Float64[], 0)
        end
        
        # Use and cache the result
        global_min_mz = g_min_mz
        global_max_mz = g_max_mz
        set_global_mz_range!(msi_data, global_min_mz, global_max_mz)
        println("  Global m/z range found and cached: [$(global_min_mz), $(global_max_mz)]")
        @printf "  (Pass 1 took %.2f seconds)\n" pass1_duration
    end
    
    # 2. Define Bins and precompute constants
    mz_bins = range(global_min_mz, stop=global_max_mz, length=num_bins)
    bin_step = step(mz_bins)
    inv_bin_step = 1.0 / bin_step  # Precompute reciprocal to avoid division
    min_mz = global_min_mz

    # Initialize thread-local intensity sums
    thread_local_intensity_sums = [zeros(Float64, num_bins) for _ in 1:Base.Threads.nthreads()]
    thread_local_spectra_processed = [0 for _ in 1:Base.Threads.nthreads()]

    # 3. Second Pass: Optimized binning
    pass2_start_time = time_ns()
    println("  Pass 2: Summing intensities into $num_bins bins...")
    _iterate_spectra_fast(msi_data, masked_indices === nothing ? nothing : collect(masked_indices)) do idx, mz, intensity
        thread_id = Base.Threads.threadid()
        thread_local_spectra_processed[thread_id] += 1
        if isempty(mz)
            return
        end
        
        # Use views to avoid bounds checks on the input arrays
        mz_view = mz
        intensity_view = intensity
        
        # Manual binning without clamp/round for SIMD
        @inbounds @simd for i in eachindex(mz_view)
            # Calculate raw bin index (much faster than round+clamp)
            raw_index = (mz_view[i] - min_mz) * inv_bin_step + 1.0
            bin_index = trunc(Int, raw_index)
            
            # Branchless version using clamp
            final_index = clamp(bin_index, 1, num_bins)
            thread_local_intensity_sums[thread_id][final_index] += intensity_view[i]
        end
    end

    # Combine thread-local results
    intensity_sum = zeros(Float64, num_bins)
    num_spectra_processed = 0
    for i in 1:Base.Threads.nthreads()
        intensity_sum .+= thread_local_intensity_sums[i]
        num_spectra_processed += thread_local_spectra_processed[i]
    end
    pass2_duration = (time_ns() - pass2_start_time) / 1e9

    total_duration = (time_ns() - total_start_time) / 1e9
    println("\n--- Profiling (Post-Optimization) ---")
    @printf "  Pass 2 (I/O + Binning):   %.2f seconds\n" pass2_duration
    @printf "  Total Function Time:      %.2f seconds\n" total_duration
    println("----------------------------------------\n")

    println("Total spectrum calculation complete.")
    return (collect(mz_bins), intensity_sum, num_spectra_processed)
end

"""
    get_total_spectrum(msi_data::MSIData; num_bins::Int=2000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the sum of all spectra in the dataset by binning.
This function dispatches to a specialized implementation based on the file type
(.imzML or .mzML) for optimal performance.

Returns a tuple containing two vectors: the binned m/z axis and the summed intensities.

# Arguments

- `msi_data::MSIData`: The MSIData object.
- `num_bins::Int`: The number of bins to use for the total spectrum.
- `mask_path::Union{String, Nothing}`: The path to the mask file.

# Returns
- `(mz_range, total_spectrum)`: A tuple containing the m/z range and the total spectrum.

# Example

```julia
get_total_spectrum(msi_data)
```
"""
function get_total_spectrum(msi_data::MSIData; num_bins::Int=2000, mask_path::Union{String, Nothing}=nothing)::Tuple{Vector{Float64}, Vector{Float64}, Int}
    local masked_indices::Union{Set{Int}, Nothing} = nothing
    if mask_path !== nothing
        mask_matrix = load_and_prepare_mask(mask_path, msi_data.image_dims)
        masked_indices = get_masked_spectrum_indices(msi_data, mask_matrix)
        println("Calculating total spectrum for masked region (mask from: $(mask_path))")
    end

    return get_total_spectrum_core(msi_data, num_bins=num_bins, masked_indices=masked_indices)
end

"""
    get_average_spectrum(msi_data::MSIData; num_bins::Int=2000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the average of all spectra in the dataset by binning.
This is effectively the total spectrum divided by the number of spectra.

Returns a tuple containing two vectors: the binned m/z axis and the averaged intensities.

# Arguments

- `msi_data::MSIData`: The MSIData object.
- `num_bins::Int`: The number of bins to use for the total spectrum.
- `mask_path::Union{String, Nothing}`: The path to the mask file.

# Returns
- `(mz_range, average_spectrum)`: A tuple containing the m/z range and the average spectrum.

# Example

```julia
get_average_spectrum(msi_data)
```
"""
function get_average_spectrum(msi_data::MSIData; num_bins::Int=2000, mask_path::Union{String, Nothing}=nothing)::Tuple{Vector{Float64}, Vector{Float64}}
    mz_bins, intensity_sum, num_spectra_processed = get_total_spectrum(msi_data, num_bins=num_bins, mask_path=mask_path)

    if isempty(intensity_sum) || num_spectra_processed == 0
        @warn "Cannot calculate average spectrum: no spectra were processed."
        return (mz_bins, Float64[])
    end

    average_intensity = intensity_sum ./ num_spectra_processed
    
    return (mz_bins, average_intensity)
end

# --- Iterator Implementation --- #

"""
    MSIDataIterator

A stateful iterator for sequentially accessing spectra from an `MSIData` object.
This is created by the `IterateSpectra` function.

# Arguments

- `data::MSIData`: The MSIData object.

# Returns
- `MSIDataIterator`: An iterator that yields each spectrum.

# Example

```julia
IterateSpectra(data)
```
"""
struct MSIDataIterator
    data::MSIData
end

"""
    IterateSpectra(data::MSIData)

Returns an iterator that yields each spectrum, processing the file sequentially.
This iterator is useful for processing all spectra in a loop and benefits from
the caching implemented in `GetSpectrum`.

This function is the core of the "Event-driven" access pattern.

# Arguments

- `data::MSIData`: The MSIData object.

# Returns
- `MSIDataIterator`: An iterator that yields each spectrum.

# Example

```julia
IterateSpectra(data)
```
"""
function IterateSpectra(data::MSIData)
    return MSIDataIterator(data)
end

# Define the iteration interface for the custom iterator
Base.length(it::MSIDataIterator) = length(it.data.spectra_metadata)
Base.eltype(it::MSIDataIterator) = Tuple{Int, Tuple{Vector, Vector}} # Index, (mz, intensity)

"""
    Base.iterate(it::MSIDataIterator, state=1)

Advances the iterator, returning the next spectrum in the sequence. It calls
`GetSpectrum`, which means the iteration process is cached according to the
`MSIData` object's cache settings.

# Arguments
- `it`: The `MSIDataIterator` instance.
- `state`: The index of the next spectrum to retrieve.

# Returns
- A tuple `((index, spectrum), next_state)` or `nothing` if iteration is complete.
"""
function Base.iterate(it::MSIDataIterator, state=1)
    if state > length(it.data.spectra_metadata)
        return nothing # End of iteration
    end

    # GetSpectrum uses the cache, so iteration is also cached
    spectrum = GetSpectrum(it.data, state)
    
    # Yield the index and the spectrum
    return ((state, spectrum), state + 1)
end


# --- High-performance Internal Iterator --- #

"""
    read_compressed_array(data::MSIData, io::IO, asset::SpectrumAsset, ::Type{T}) where {T}

Reads a single data array (m/z or intensity) from an `.ibd` file stream,
handling both compressed and uncompressed data.

This is an internal function designed for high-performance reading and decompressing of binary arrays.
Uses type parameters and buffer pooling to minimize allocations and maximize speed.

# Arguments
- `data`: The `MSIData` object.
- `io`: The IO stream of the `.ibd` file.
- `asset`: The `SpectrumAsset` for the array.
- `::Type{T}`: The target format of the data.

# Returns
- A `Vector{T}` containing the data.

# Throws
- An error if zlib decompression fails.
"""
function read_compressed_array(data::MSIData, io::IO, asset::SpectrumAsset, ::Type{T}) where {T}
    # Add validation before seeking
    if asset.offset < 0 || asset.offset >= filesize(io)
        throw(FileFormatError("Invalid asset offset: $(asset.offset) for file size $(filesize(io))"))
    end
    # Optimization: Use Mmap if available to avoid seek and copy
    mmap_data = (data.source isa ImzMLSource) ? data.source.mmap_data : nothing

    if asset.is_compressed
        local decompressed_view
        if mmap_data !== nothing
            # Zero-copy access to the compressed segment
            # Base64 should be read using a view as well
            compressed_view = view(mmap_data, (asset.offset + 1):(asset.offset + asset.encoded_length))
            decompressed_view = Libz.inflate(compressed_view)
        else
            # Fallback to standard IO
            seek(io, asset.offset)
            compressed_bytes_buffer = get_buffer!(data.buffer_pool, Int(asset.encoded_length))
            try
                readbytes!(io, compressed_bytes_buffer, asset.encoded_length)
                decompressed_view = Libz.inflate(view(compressed_bytes_buffer, 1:asset.encoded_length))
            finally
                release_buffer!(data.buffer_pool, compressed_bytes_buffer)
            end
        end
        
        local array
        try
            # Pre-allocate the typed output array
            n_elements = length(decompressed_view) ÷ sizeof(T)
            array = Vector{T}(undef, n_elements)
            
            # Use unsafe_copyto! for zero-overhead copy into the typed array
            unsafe_copyto!(reinterpret(Ptr{UInt8}, pointer(array)), pointer(decompressed_view), length(decompressed_view))
            
        catch e
            @error "ZLIB DECOMPRESSION FAILED at offset $(asset.offset)"
            rethrow(e)
        end
        
        return array
    else
        # Read uncompressed data
        if mmap_data !== nothing
            # SAFETY: Check address alignment (sizeof(T) must divide asset.offset)
            # Since mmap_data itself is page-aligned, we only check the offset.
            alignment = sizeof(T)
            if asset.offset % alignment == 0
                # ZERO-COPY Path
                end_pos = asset.offset + asset.encoded_length * alignment
                raw_view = view(mmap_data, (asset.offset + 1):end_pos)
                return Vector{T}(reinterpret(T, raw_view))
            else
                # ALIGNMENT FALLBACK: Memory-to-memory copy (safer than reinterpret)
                array = Vector{T}(undef, asset.encoded_length)
                # Raw copy from mmap to vector
                n_bytes = asset.encoded_length * alignment
                unsafe_copyto!(reinterpret(Ptr{UInt8}, pointer(array)), pointer(mmap_data, asset.offset + 1), n_bytes)
                return array
            end
        else
            # Fallback to standard IO
            seek(io, asset.offset)
            array = Vector{T}(undef, asset.encoded_length)
            read!(io, array)
            return array
        end
    end
end


"""
    _iterate_uncompressed_fast(f::Function, data::MSIData, source::ImzMLSource)

A highly optimized iterator for uncompressed `.imzML` data.

It pre-allocates large buffers for m/z and intensity arrays and reuses them
for each spectrum by creating views. This minimizes memory allocations and
is significantly faster for bulk processing than reading spectra one by one.

# Arguments
- `f`: A function to execute for each spectrum, with the signature `f(index, mz_view, int_view)`.
- `data`: The `MSIData` object.
- `source`: The `ImzMLSource`.

# Returns
- `nothing`

# Example

```julia
_iterate_uncompressed_fast(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""

"""
    _iterate_compressed_fast(f::Function, data::MSIData, source::ImzMLSource)

An iterator for `.imzML` datasets that contain compressed spectra.

This function iterates through each spectrum and reads its data individually,
decompressing it on the fly if necessary. Because the size of decompressed
data is not known in advance, this path cannot use the buffer-reuse optimization
and will be slower and allocate more memory than `_iterate_uncompressed_fast`.

# Arguments
- `f`: A function to execute for each spectrum, with the signature `f(index, mz_array, int_array)`.
- `data`: The `MSIData` object.
- `source`: The `ImzMLSource`.

# Returns
- `nothing`

# Example

```julia
_iterate_compressed_fast(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function _iterate_compressed_fast(f::Function, data::MSIData, source::ImzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    # Path for datasets containing at least one compressed spectrum.
    # Optimized to minimize allocations by reusing buffers.
    
    # Determine which indices to iterate over
    spectrum_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate

    # Pre-allocate large enough buffers for the expected maximum number of points
    # We estimate based on metadata if possible, or grow dynamically
    max_encoded = maximum(meta -> max(meta.mz_asset.encoded_length, meta.int_asset.encoded_length), data.spectra_metadata)
    
    # Heuristic: decompressed size is usually larger. We'll start with 10x and grow if needed.
    # But read_compressed_array currently returns a Vector, so we'll need to modify it 
    # to accept an optional target buffer.
    
    for i in spectrum_indices
        meta = data.spectra_metadata[i]
        
        if meta.mz_asset.encoded_length == 0 && meta.int_asset.encoded_length == 0
            f(i, source.mz_format[], source.intensity_format[])
            continue
        end

        # For compressed data, we currently allocate new arrays per spectrum.
        # To truly minimize allocations, we'd need read_compressed_array! (in-place version).
        # For now, let's ensure we are at least using the optimized type-stable version.
        handle = get_handle(source)
        mz_array = read_compressed_array(data, handle, meta.mz_asset, source.mz_format)
        intensity_array = read_compressed_array(data, handle, meta.int_asset, source.intensity_format)
        
        mz_array .= ltoh.(mz_array)
        intensity_array .= ltoh.(intensity_array)
        
        f(i, mz_array, intensity_array)
    end
end

"""
    _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource)

Internal implementation of the fast iterator for `.imzML` files. It reads
data directly from the `.ibd` file stream.

This function acts as a dispatcher:
- If any spectrum is compressed, it uses a slower path that decompresses each spectrum individually.
- If all spectra are uncompressed, it uses a highly optimized path that reuses pre-allocated 
  buffers to minimize memory allocations and overhead.

# Arguments

- `f::Function`: A function to execute for each spectrum, with the signature `f(index, mz_array, int_array)`.
- `data::MSIData`: The `MSIData` object.
- `source::ImzMLSource`: The `ImzMLSource`.

# Returns
- `nothing`

# Example

```julia
_iterate_spectra_fast_impl(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    if isempty(data.spectra_metadata)
        return
    end
    
    # Use cached compression status for zero-allocation dispatch
    if source.is_any_compressed
        _iterate_compressed_fast(f, data, source, indices_to_iterate)
    else
        _iterate_uncompressed_fast(f, data, source, indices_to_iterate)
    end
end

"""
    _iterate_spectra_fast_impl(f::Function, data::MSIData, source::MzMLSource)

Internal implementation of the fast iterator for `.mzML` files. It iterates
through each spectrum, decodes the Base64 data on the fly, and calls the
provided function. It bypasses the `GetSpectrum` cache to avoid storing
all decoded spectra in memory.

# Arguments

- `f::Function`: A function to execute for each spectrum, with the signature `f(index, mz_array, int_array)`.
- `data::MSIData`: The `MSIData` object.
- `source::MzMLSource`: The `MzMLSource`.

# Returns
- `nothing`

# Example

```julia
_iterate_spectra_fast_impl(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::MzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    spectrum_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate
    indices_with_offsets = [(i, data.spectra_metadata[i].mz_asset.offset) for i in spectrum_indices]
    sort!(indices_with_offsets, by = x -> x[2])

    handle = get_handle(source)

    for (i, _) in indices_with_offsets
        meta = data.spectra_metadata[i]
        # For MzML, we still have some allocations due to Base64 decoding, 
        # but we use the thread-local handle.
        mz = read_binary_vector(data, handle, meta.mz_asset)
        intensity = read_binary_vector(data, handle, meta.int_asset)
        f(i, mz, intensity)
    end
end

function _iterate_uncompressed_fast(f::Function, data::MSIData, source::ImzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    spectrum_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate
    
    # Intialize local buffers for this thread's sequential iteration
    mz_buf = Vector{Float64}()
    int_buf = Vector{Float64}()

    for i in spectrum_indices
        meta = data.spectra_metadata[i]
        
        # Use Mmap views if available (zero-copy if LE)
        if source.mmap_data !== nothing
            # Optimized Mmap path (same as read_spectrum_from_disk but potentially avoiding copies)
            mz, intensity = read_spectrum_from_disk(source, meta)
            f(i, mz, intensity)
        else
            # Read into our loop buffers to avoid continuous allocation
            handle = get_handle(source)
            
            # Resize buffers if necessary (minimal reallocation)
            resize!(mz_buf, meta.mz_asset.encoded_length)
            resize!(int_buf, meta.int_asset.encoded_length)
            
            seek(handle, meta.mz_asset.offset)
            read!(handle, mz_buf)
            seek(handle, meta.int_asset.offset)
            read!(handle, int_buf)
            
            # Convert in-place if possible
            mz_buf .= ltoh.(mz_buf)
            int_buf .= ltoh.(int_buf)
            
            f(i, mz_buf, int_buf)
        end
    end
end

"""
    _iterate_spectra_fast_serial(f::Function, data::MSIData, indices_to_iterate=nothing)

An internal, high-performance iterator for bulk processing that bypasses the `GetSpectrum` cache.
It provides direct, sequential access to spectral data, making it the most efficient
method for operations that need to process many or all spectra (e.g., calculating
a total spectrum, pre-computing analytics, or exporting data).

This function dispatches to a specialized implementation based on the data source type
(`ImzMLSource` or `MzMLSource`) and data characteristics (e.g., compressed vs.
uncompressed) to maximize performance.

# Arguments
- `f`: A function to call for each spectrum, with the signature `f(index, mz, intensity)`.
- `data`: The `MSIData` object.
- `indices_to_iterate`: An optional `AbstractVector{Int}` specifying a subset of
  spectrum indices to iterate over. If `nothing`, it iterates over all spectra.

# Warning
This is a low-level function that bypasses the public `GetSpectrum` API. It does not
use the cache and is not guaranteed to be thread-safe if the same `MSIData` object
is accessed from multiple threads concurrently.
"""
function _iterate_spectra_fast_serial(f::Function, data::MSIData, indices_to_iterate::Union{AbstractVector{Int}, Nothing}=nothing)
    # Dispatch to the correct implementation based on the source type
    _iterate_spectra_fast_impl(f, data, data.source, indices_to_iterate)
end

"""
    _iterate_spectra_fast(f::Function, data::MSIData, indices_to_iterate=nothing)

A smart dispatcher for internal, high-performance iteration over spectra.
It automatically chooses between serial and parallel execution based on the
number of available threads (`Base.Threads.nthreads()`).

# Arguments
- `f`: A function to call for each spectrum, with the signature `f(index, mz, intensity)`.
- `data`: The `MSIData` object.
- `indices_to_iterate`: An optional `AbstractVector{Int}` specifying a subset of
  spectrum indices to iterate over. If `nothing`, it iterates over all spectra.

# Returns
- `nothing`

# Example

```julia
_iterate_spectra_fast(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function _iterate_spectra_fast(f::Function, data::MSIData, indices_to_iterate::Union{AbstractVector{Int}, Nothing}=nothing)
    # Default to all indices if not specified
    all_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate
    
    if Base.Threads.nthreads() > 1 && length(all_indices) > 100 # Heuristic: only parallelize for enough work
        _iterate_spectra_fast_parallel(f, data, all_indices)
    else
        _iterate_spectra_fast_serial(f, data, all_indices)
    end
end

"""
    _iterate_spectra_fast_threadsafe(f::Function, data::MSIData, indices::AbstractVector{Int})

Thread-safe version of the fast iterator that uses the cache lock to ensure exclusive file access.
This is safe to use when multiple threads might access the same MSIData object concurrently.

# Arguments
- `f`: Function to call for each spectrum with signature `f(index, mz, intensity)`
- `data`: The MSIData object
- `indices`: Specific spectrum indices to iterate over

# Performance
- Slower than non-thread-safe version due to locking
- Safe for concurrent access to same MSIData object
- Good for small to medium datasets where parallel overhead isn't justified
"""
function _iterate_spectra_fast_threadsafe(f::Function, data::MSIData, indices::AbstractVector{Int})
    # Use the cache lock to ensure exclusive access to file
    lock(data.cache_lock) do
        _iterate_spectra_fast(f, data, indices)
    end
end

"""
    _iterate_spectra_fast_parallel(f::Function, data::MSIData, indices::AbstractVector{Int})

Parallel version that creates thread-local file handles for maximum performance.
Each thread gets its own file handle, eliminating contention.

# Arguments
- `f`: Function to call for each spectrum with signature `f(index, mz, intensity)`
- `data`: The MSIData object  
- `indices`: Specific spectrum indices to iterate over

# Performance
- Fastest for large datasets on multi-core systems
- Eliminates file I/O contention between threads
- Higher memory usage due to multiple file handles
- Best for bulk processing operations
"""
function _iterate_spectra_fast_parallel(f::Function, data::MSIData, indices::AbstractVector)
    n_total = length(indices)
    n_threads = Base.Threads.nthreads()
    
    # Manual chunking to avoid allocations of Iterators.partition and collect
    Base.Threads.@threads for t in 1:n_threads
        start_idx = ((t - 1) * n_total ÷ n_threads) + 1
        end_idx = (t * n_total) ÷ n_threads
        
        if start_idx <= end_idx
            chunk = view(indices, start_idx:end_idx)
            _iterate_spectra_fast_impl(f, data, data.source, chunk)
        end
    end
end

"""
    _iterate_spectra_fast_parallel(f::Function, data::MSIData)

Parallel version that processes all spectra in the dataset.
Convenience wrapper for the indexed version.

# Arguments

- `f::Function`: A function to execute for each spectrum, with the signature `f(index, mz_array, int_array)`.
- `data::MSIData`: The `MSIData` object.

# Returns
- `nothing`

# Example

```julia
_iterate_spectra_fast_parallel(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function _iterate_spectra_fast_parallel(f::Function, data::MSIData)
    all_indices = collect(1:length(data.spectra_metadata))
    _iterate_spectra_fast_parallel(f, data, all_indices)
end

"""
    _iterate_spectra_fast_threadsafe(f::Function, data::MSIData)

Thread-safe version that processes all spectra in the dataset.
Convenience wrapper for the indexed version.

# Arguments

- `f::Function`: A function to execute for each spectrum, with the signature `f(index, mz_array, int_array)`.
- `data::MSIData`: The `MSIData` object.

# Returns
- `nothing`

# Example

```julia
_iterate_spectra_fast_threadsafe(data, 1) do mz, intensity
    # Process the spectrum here
end
```
"""
function _iterate_spectra_fast_threadsafe(f::Function, data::MSIData)
    all_indices = collect(1:length(data.spectra_metadata))
    _iterate_spectra_fast_threadsafe(f, data, all_indices)
end

"""
    get_bloom_filter(data::MSIData, index::Int)::BloomFilter{Int}

Retrieves the Bloom filter for a specific spectrum by its index.
If the Bloom filter has not yet been created (due to lazy initialization),
it will be created on-demand and cached.

# Arguments
- `data`: The `MSIData` object.
- `index`: The linear index of the spectrum.

# Returns
- A `BloomFilter{Int}` for the specified spectrum.
"""
function get_bloom_filter(data::MSIData, index::Int)::BloomFilter{Int}
    # Double-checked locking for resilience if not precomputed
    if data.bloom_filters === nothing || data.bloom_filters[index] === nothing
        lock(data.cache_lock) do
            if data.bloom_filters === nothing
                data.bloom_filters = Vector{Union{BloomFilter{Int}, Nothing}}(undef, length(data.spectra_metadata))
                fill!(data.bloom_filters, nothing)
            end
            
            if data.bloom_filters[index] === nothing
                mz, _ = GetSpectrum(data, index)
                data.bloom_filters[index] = create_bloom_filter_for_spectrum(mz)
            end
        end
    end
    return data.bloom_filters[index]
end
