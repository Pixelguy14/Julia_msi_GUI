# src/MSIData.jl

using Base64, Libz, Serialization, Printf, DataFrames, Base.Threads, StatsBase
using Mmap

const FILE_HANDLE_LOCK = ReentrantLock()

"""
    calculate_optimal_analytics_chunk_size(total_spectra::Int, system_caps::NamedTuple) -> Int

Calculates an optimal chunk size for analytics processing based on total number of spectra
and detected system capabilities (total memory, number of CPU threads).

The heuristic aims to balance memory usage and parallelism:
- Smaller chunks for less memory or fewer threads.
- Larger chunks for more memory or more threads, but with a practical upper limit.
"""
function calculate_optimal_analytics_chunk_size(total_spectra::Int, system_caps::NamedTuple)::Int
    # Base chunk size - can be refined with more profiling
    base_chunk_size = 5000

    # Adjust based on available RAM (heuristic: more RAM, larger chunks)
    # Assume 1GB RAM can handle a certain chunk size comfortably
    # E.g., for 32GB RAM, allow 2x more than for 16GB
    memory_factor = max(1.0, system_caps.total_memory_gb / 16.0) # Scale based on 16GB baseline

    # Adjust based on CPU threads (heuristic: more threads, potentially more, smaller chunks to keep threads busy, or fewer large chunks)
    # A simple approach for now: scale with sqrt of threads to avoid overly large chunks with many threads
    thread_factor = max(1.0, sqrt(system_caps.num_cpu_threads / 4.0)) # Scale based on 4-thread baseline

    # Combine factors, with a practical minimum and maximum
    dynamic_chunk_size = round(Int, base_chunk_size * memory_factor * thread_factor)
    
    # Ensure chunk size is within reasonable bounds
    min_chunk_size = 500
    max_chunk_size = 20000

    return clamp(dynamic_chunk_size, min_chunk_size, max_chunk_size)

end

"""
    ThreadSafeFileHandle

Wraps a file handle with thread-safe operations.
"""
mutable struct ThreadSafeFileHandle
    path::String
    handle::IO
    lock::ReentrantLock
end

function ThreadSafeFileHandle(path::String, mode::String="r")
    handle = lock(FILE_HANDLE_LOCK) do
        open(path, mode)
    end
    return ThreadSafeFileHandle(path, handle, ReentrantLock())
end

"""
    read_at!(tsfh::ThreadSafeFileHandle, a::AbstractArray, pos::Integer)

Atomically seeks to a position and reads data into an array. This is the thread-safe
way to read from a specific offset in the file.
"""
function read_at!(tsfh::ThreadSafeFileHandle, a::AbstractArray, pos::Integer)
    lock(tsfh.lock) do
        try
            seek(tsfh.handle, pos)
            read!(tsfh.handle, a)
        catch e
            error_str = sprint(showerror, e)
            if isa(e, SystemError) && occursin("Bad file descriptor", error_str)
                @warn "Bad file descriptor detected. Attempting to re-open file handle for $(tsfh.path) and retry."
                close(tsfh.handle) # Close the problematic handle
                tsfh.handle = open(tsfh.path, "r") # Re-open
                seek(tsfh.handle, pos) # Retry seek
                read!(tsfh.handle, a) # Retry read
            else
                rethrow(e)
            end
        end
    end
end

function Base.read(tsfh::ThreadSafeFileHandle, n::Integer)
    lock(tsfh.lock) do
        return read(tsfh.handle, n)
    end
end

function Base.close(tsfh::ThreadSafeFileHandle)
    close(tsfh.handle)
end

function Base.isopen(tsfh::ThreadSafeFileHandle)
    isopen(tsfh.handle)
end

function Base.position(tsfh::ThreadSafeFileHandle)
    lock(tsfh.lock) do
        return position(tsfh.handle)
    end
end

function Base.filesize(tsfh::ThreadSafeFileHandle)
    lock(tsfh.lock) do
        return filesize(tsfh.handle)
    end
end

# Overload Base.seek for ThreadSafeFileHandle
function Base.seek(tsfh::ThreadSafeFileHandle, pos::Integer)
    lock(tsfh.lock) do
        seek(tsfh.handle, pos)
    end
end

# Overload Base.read! for ThreadSafeFileHandle
function Base.read!(tsfh::ThreadSafeFileHandle, a::AbstractArray)
    lock(tsfh.lock) do
        read!(tsfh.handle, a)
    end
end

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
    ibd_handle::Union{IO, ThreadSafeFileHandle, Vector{UInt8}}
    mz_format::Type
    intensity_format::Type
end

function ImzMLSource(ibd_path::String, mz_format::Type, intensity_format::Type, use_mmap::Bool)
    if use_mmap
        println("Using memory-mapping for .ibd file.")
        try
            mmap_array = Mmap.mmap(ibd_path, Vector{UInt8})
            return ImzMLSource(mmap_array, mz_format, intensity_format)
        catch e
            @warn "Memory mapping failed, falling back to file handles: $e"
        end
    end
    
    println("Using file handles for .ibd file.")
    handle = ThreadSafeFileHandle(ibd_path, "r")
    return ImzMLSource(handle, mz_format, intensity_format)
end

"""
    MzMLSource <: MSDataSource

A data source for `.mzML` files, holding a handle to the `.mzML` file itself
(which contains the binary data encoded in Base64) and the expected data formats.
"""
struct MzMLSource <: MSDataSource
    file_handle::Union{IO, ThreadSafeFileHandle}
    mz_format::Type
    intensity_format::Type
end

"""
    SpectrumAsset

Contains metadata for a single binary data array (m/z or intensity) within a spectrum.

# Fields
- `format`: The data type of the elements (e.g., `Float32`, `Int64`).
- `is_compressed`: A boolean flag indicating if the data is compressed (e.g., with zlib).
- `offset`: The byte offset of the data within the file (`.ibd` for imzML, `.mzML` for mzML).
- `encoded_length`:
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
- `type`: e.g., :sample, :blank, :qc
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
    MSIData

The primary object for interacting with mass spectrometry data. It provides a unified
interface for both `.mzML` and `.imzML` files and includes an LRU cache for
efficient repeated access to spectra.

# Fields
- `source`: The underlying `MSDataSource` (`ImzMLSource` or `MzMLSource`).
- `spectra_metadata`: A vector of `SpectrumMetadata` for all spectra in the file.
- `instrument_metadata`: Metadata about the instrument configuration.
- `image_dims`: A tuple `(width, height)` of the spatial dimensions (for imzML).
- `coordinate_map`: A matrix mapping `(x, y)` coordinates to a linear spectrum index (for imzML).
- `cache`: A dictionary holding cached spectra, mapping index to `(mz, intensity)`.
- `cache_order`: A vector of indices tracking usage for the LRU cache policy.
- `cache_size`: The maximum number of spectra to store in the cache.
- `cache_lock`: A `ReentrantLock` to ensure thread-safe access to the cache.
- `global_min_mz`: Cached global minimum m/z value across all spectra.
- `global_max_mz`: Cached global maximum m/z value across all spectra.
- `spectrum_stats_df`: A `DataFrame` containing pre-computed per-spectrum analytics (e.g., TIC, BPI).
- `bloom_filters`: Lazy-initialized Bloom filters for fast spectrum content checking.
- `analytics_ready`: Flag indicating if analytics have been pre-computed.
- `preprocessing_hints`: Dictionary for auto-determined parameters.
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
    buffer_pool::SimpleBufferPool

    # Pre-computed analytics/metadata - use Threads.Atomic for compatibility
    global_min_mz::Threads.Atomic{Float64}
    global_max_mz::Threads.Atomic{Float64}
    spectrum_stats_df::Union{DataFrame, Nothing}
    bloom_filters::Union{Vector{<:BloomFilter}, Nothing}
    analytics_ready::AtomicFlag
    preprocessing_hints::Union{Dict{Symbol, Any}, Nothing} # For auto-determined parameters

    function MSIData(source, metadata, instrument_meta, dims, coordinate_map, cache_size)
        obj = new(source, metadata, instrument_meta, dims, coordinate_map, 
                  Dict(), [], cache_size, ReentrantLock(),
                  SimpleBufferPool(),
                  Threads.Atomic{Float64}(Inf), Threads.Atomic{Float64}(-Inf), 
                  nothing, nothing, AtomicFlag(), nothing)
        
        # Ensure file handles are closed when the object is garbage collected
        finalizer(obj) do o
            if o.source isa ImzMLSource && isopen(o.source.ibd_handle)
                close(o.source.ibd_handle)
            elseif o.source isa MzMLSource && isopen(o.source.file_handle)
                close(o.source.file_handle)
            end
        end
        return obj
    end
end

# Thread-safe state management
function set_global_mz_range!(data::MSIData, min_mz::Float64, max_mz::Float64)
    Threads.atomic_xchg!(data.global_min_mz, min_mz)
    Threads.atomic_xchg!(data.global_max_mz, max_mz)
end

function get_global_mz_range(data::MSIData)
    # Atomic read using atomic_add! with zero
    min_mz = Threads.atomic_add!(data.global_min_mz, 0.0)
    max_mz = Threads.atomic_add!(data.global_max_mz, 0.0)
    return (min_mz, max_mz)
end

function set_analytics_data!(data::MSIData, stats_df::DataFrame, bloom_filters::Vector)
    lock(data.cache_lock) do
        data.spectrum_stats_df = stats_df
        data.bloom_filters = bloom_filters
    end
end

function get_bloom_filters(data::MSIData)
    lock(data.cache_lock) do
        return data.bloom_filters
    end
end

function get_spectrum_stats(data::MSIData)
    lock(data.cache_lock) do
        return data.spectrum_stats_df
    end
end

"""
    Base.close(data::MSIData)

Explicitly closes the file handles associated with the `MSIData` object and clears the spectrum cache.
It is good practice to call this method when you are finished with an `MSIData` object to release resources immediately.
"""
function Base.close(data::MSIData)
    if data.source isa ImzMLSource && isopen(data.source.ibd_handle)
        close(data.source.ibd_handle)
    elseif data.source isa MzMLSource && isopen(data.source.file_handle)
        close(data.source.file_handle)
    end

    # Clear cache
    empty!(data.cache)
    empty!(data.cache_order)

    println("MSIData object closed and resources released.")
end

"""
    warm_cache(data::MSIData, indices::AbstractVector{Int})

Pre-loads a specified list of spectra into the cache. This is useful for warming up
the cache with frequently accessed spectra to improve performance for subsequent
interactive analysis.
"""
function warm_cache(data::MSIData, indices::AbstractVector{Int})
    println("Warming cache with $(length(indices)) spectra...")
    for idx in indices
        # Calling GetSpectrum will load the data and place it in the cache
        GetSpectrum(data, idx)
    end
    println("Cache warming complete.")
end

"""
    get_memory_usage(data::MSIData)

Calculates and returns a summary of the memory currently used by the MSIData object,
including the spectrum cache and metadata.

# Returns
- A `NamedTuple` with memory usage in bytes for different components.
"""
function get_memory_usage(data::MSIData)
    cache_bytes = 0
    lock(data.cache_lock) do
        # This is an approximation; `sizeof` on a dictionary is not fully representative,
        # but it's much better to sum the size of the actual vectors.
        for (mz, intensity) in values(data.cache)
            cache_bytes += sizeof(mz) + sizeof(intensity)
        end
    end

    metadata_bytes = sizeof(data.spectra_metadata)
    coordinate_map_bytes = data.coordinate_map === nothing ? 0 : sizeof(data.coordinate_map)

    total_bytes = cache_bytes + metadata_bytes + coordinate_map_bytes

    return (
        total_mb = total_bytes / 1024^2,
        cache_mb = cache_bytes / 1024^2,
        metadata_mb = metadata_bytes / 1024^2,
        coordinate_map_mb = coordinate_map_bytes / 1024^2
    )
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
    if asset.offset < 0 || asset.offset >= filesize(io)
        throw(FileFormatError("Invalid asset offset: $(asset.offset) for file size $(filesize(io))"))
    end
    
    seek(io, asset.offset)
    raw_b64 = read(io, asset.encoded_length)
    
    # Use String directly to avoid intermediate allocations
    b64_string = String(raw_b64)
    
    local decoded_bytes::Vector{UInt8}
    
    if asset.is_compressed
        # Direct Base64 decode to temporary, then decompress
        temp_decoded = Base64.base64decode(b64_string)
        decoded_bytes = Libz.inflate(temp_decoded)
    else
        # Direct Base64 decode
        decoded_bytes = Base64.base64decode(b64_string)
    end
    
    # Calculate number of elements
    n_elements = length(decoded_bytes) ÷ sizeof(asset.format)

    if n_elements * sizeof(asset.format) != length(decoded_bytes)
        throw(FileFormatError("Size of decoded byte array is not a multiple of the element size."))
    end
    
    # Reinterpret the byte array as an array of the target type. This does not copy.
    reinterpreted_array = reinterpret(asset.format, decoded_bytes)
    
    # Allocate the final output array and convert byte order while copying.
    out_array = [ltoh(x) for x in reinterpreted_array]
    
    return out_array
end

function read_binary_vector(data::MSIData, ts_handle::ThreadSafeFileHandle, asset::SpectrumAsset)
    lock(ts_handle.lock) do
        return read_binary_vector(data, ts_handle.handle, asset)
    end
end

# Overload for different source types
"""
    read_uncompressed_array(mmap_array::Vector{UInt8}, asset::SpectrumAsset, format::Type)

Reads a single spectrum's m/z and intensity arrays directly from the `.ibd`
binary file for an `.imzML` dataset. It handles reading the raw binary data
and converting it from little-endian to the host's native byte order.

# Arguments
- `source`: An `ImzMLSource` containing the file handle and data formats.
- `meta`: The `SpectrumMetadata` for the spectrum to be read.

# Returns
- A tuple `(mz, intensity)` containing the two requested data arrays.
"""
function read_uncompressed_array(mmap_array::Vector{UInt8}, asset::SpectrumAsset, format::Type)
    # For uncompressed imzML, encoded_length is number of elements
    n_elements = asset.encoded_length
    byte_length = n_elements * sizeof(format)
    
    # Check bounds
    if asset.offset + byte_length > length(mmap_array)
        throw(BoundsError(mmap_array, asset.offset+1:asset.offset+byte_length))
    end
    
    # Create a view into the mmapped array
    byte_view = @view mmap_array[asset.offset+1:asset.offset+byte_length]
    
    # Reinterpret as the desired type
    reinterpreted = reinterpret(format, byte_view)
    
    # Convert endianness (imzML is little-endian)
    return [ltoh(x) for x in reinterpreted]
end

function read_spectrum_from_disk(source::ImzMLSource, meta::SpectrumMetadata)
    if source.ibd_handle isa Vector{UInt8}
        # Fast path: read from mmapped array
        mz = read_uncompressed_array(source.ibd_handle, meta.mz_asset, source.mz_format)
        intensity = read_uncompressed_array(source.ibd_handle, meta.int_asset, source.intensity_format)
        return mz, intensity
    else
        # For imzML, the binary data is raw, not base64 encoded.
        # The `encoded_length` field in this case holds the number of points.
        mz = Array{source.mz_format}(undef, meta.mz_asset.encoded_length)
        intensity = Array{source.intensity_format}(undef, meta.int_asset.encoded_length)

        # Use the new atomic read_at! method for thread-safety
        read_at!(source.ibd_handle, mz, meta.mz_asset.offset)
        read_at!(source.ibd_handle, intensity, meta.int_asset.offset)

        # imzML data is little-endian. Convert to host byte order.
        mz .= ltoh.(mz)
        intensity .= ltoh.(intensity)

        validate_spectrum_data(mz, intensity, meta.id) 

        return mz, intensity
    end # This ends the else block
end # This ends the function

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
    # For mzML, data is Base64 encoded within the XML
    mz = read_binary_vector(data, source.file_handle, meta.mz_asset)
    intensity = read_binary_vector(data, source.file_handle, meta.int_asset)
    
    mz_f64, intensity_f64 = Float64.(mz), Float64.(intensity)

    validate_spectrum_data(mz, intensity, meta.id) 

    return mz, intensity
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
        thread_local_min_mz = [Inf for _ in 1:Threads.nthreads()]
        thread_local_max_mz = [-Inf for _ in 1:Threads.nthreads()]

        # Initialize vectors for per-spectrum stats
        tics = Vector{Float64}(undef, num_spectra)
        bpis = Vector{Float64}(undef, num_spectra)
        bp_mzs = Vector{Float64}(undef, num_spectra)
        num_points = Vector{Int}(undef, num_spectra)
        min_mzs = Vector{Float64}(undef, num_spectra)
        max_mzs = Vector{Float64}(undef, num_spectra)
        modes = Vector{SpectrumMode}(undef, num_spectra)
        
        # Don't precompute Bloom filters by default - create them lazily
        bloom_filters = nothing

        # Process in chunks to reduce memory pressure
        system_caps = detect_system_capabilities()
        chunk_size = calculate_optimal_analytics_chunk_size(num_spectra, system_caps)
        
        for chunk_start in 1:chunk_size:num_spectra
            chunk_end = min(chunk_start + chunk_size - 1, num_spectra)
            chunk_range = chunk_start:chunk_end
            
            println("Processing chunk $chunk_start - $chunk_end / $num_spectra")
            
            # Process current chunk
            _iterate_spectra_fast(msi_data, collect(chunk_range)) do idx, mz, intensity
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
                thread_id = Threads.threadid()
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
            Threads.atomic_xchg!(msi_data.global_min_mz, g_min_mz)
            Threads.atomic_xchg!(msi_data.global_max_mz, g_max_mz)
        else
            # Fallback: calculate from first spectrum
            try
                if num_spectra > 0
                    mz, _ = GetSpectrum(msi_data, 1)
                    if !isempty(mz)
                        Threads.atomic_xchg!(msi_data.global_min_mz, minimum(mz))
                        Threads.atomic_xchg!(msi_data.global_max_mz, maximum(mz))
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
    get_total_spectrum_imzml(msi_data::MSIData; num_bins::Int=2000)

Internal function to calculate the total spectrum for an `.imzML` file.

It uses a fast, two-pass approach:
1. The first pass finds the global m/z range across all spectra.
2. The second pass sums intensities into a pre-defined number of bins.

This function is highly optimized for `.imzML` by leveraging direct binary
reading and optimized binning logic. It is called by `get_total_spectrum`.
"""
function get_total_spectrum_imzml(msi_data::MSIData; num_bins::Int=2000, masked_indices::Union{Set{Int}, Nothing}=nothing)
    println("Calculating total spectrum for imzML (2-pass method)...")
    total_start_time = time_ns()

    local global_min_mz, global_max_mz

    if Threads.atomic_add!(msi_data.global_min_mz, 0.0) !== Inf
        println("  Using pre-computed m/z range.")
        global_min_mz = Threads.atomic_add!(msi_data.global_min_mz, 0.0)
        global_max_mz = Threads.atomic_add!(msi_data.global_max_mz, 0.0)
    else
        # 1. First Pass: Find the global m/z range by reading the fast .ibd file
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
            @warn "Could not determine a valid m/z range for imzML. All spectra might be empty."
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
    thread_local_intensity_sums = [zeros(Float64, num_bins) for _ in 1:Threads.nthreads()]
    thread_local_spectra_processed = [0 for _ in 1:Threads.nthreads()]

    # 3. Second Pass: Optimized binning
    pass2_start_time = time_ns()
    println("  Pass 2: Summing intensities into $num_bins bins...")
    _iterate_spectra_fast(msi_data, masked_indices === nothing ? nothing : collect(masked_indices)) do idx, mz, intensity
        thread_id = Threads.threadid()
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
    for i in 1:Threads.nthreads()
        intensity_sum .+= thread_local_intensity_sums[i]
        num_spectra_processed += thread_local_spectra_processed[i]
    end
    pass2_duration = (time_ns() - pass2_start_time) / 1e9

    total_duration = (time_ns() - total_start_time) / 1e9
    println("\n--- imzML Profiling (Post-Optimization) ---")
    @printf "  Pass 2 (I/O + Binning):   %.2f seconds\n" pass2_duration
    @printf "  Total Function Time:      %.2f seconds\n" total_duration
    println("----------------------------------------\n")

    println("Total spectrum calculation complete for imzML.")
    return (collect(mz_bins), intensity_sum, num_spectra_processed)
end

"""
    get_total_spectrum_mzml(msi_data::MSIData; num_bins::Int=2000)

Internal function to calculate the total spectrum for an `.mzML` file.

It uses a two-pass approach analogous to the `imzML` implementation:
1. The first pass finds the global m/z range by iterating through all spectra.
2. The second pass sums intensities into a pre-defined number of bins.

This function is called by `get_total_spectrum`.
"""
function get_total_spectrum_mzml(msi_data::MSIData; num_bins::Int=2000, masked_indices::Union{Set{Int}, Nothing}=nothing)
    println("Calculating total spectrum for mzML (2-pass method)...")
    total_start_time = time_ns()

    local global_min_mz, global_max_mz

    if Threads.atomic_add!(msi_data.global_min_mz, 0.0) !== Inf
        println("  Using pre-computed m/z range.")
        global_min_mz = Threads.atomic_add!(msi_data.global_min_mz, 0.0)
        global_max_mz = Threads.atomic_add!(msi_data.global_max_mz, 0.0)
    else
        # --- Pass 1: Find m/z range and cache it ---
        pass1_start_time = time_ns()
        println("  Pass 1: Finding global m/z range...")
        g_min_mz = Inf
        g_max_mz = -Inf
        _iterate_spectra_fast(msi_data, masked_indices === nothing ? nothing : collect(masked_indices)) do idx, mz, intensity
            if !isempty(mz)
                local_min, local_max = extrema(mz)
                g_min_mz = min(g_min_mz, local_min)
                g_max_mz = max(g_max_mz, local_max)
            end
        end
        pass1_duration = (time_ns() - pass1_start_time) / 1e9

        if !isfinite(g_min_mz)
            @warn "Could not determine a valid m/z range for mzML. All spectra might be empty."
            return (Float64[], Float64[], 0)
        end

        # Use and cache the result
        global_min_mz = g_min_mz
        global_max_mz = g_max_mz
        set_global_mz_range!(msi_data, global_min_mz, global_max_mz)

        println("  Global m/z range found and cached: [$(global_min_mz), $(global_max_mz)]")
        @printf "  (Pass 1 took %.2f seconds)\n" pass1_duration
    end

    # --- Pass 2: Bin intensities ---
    pass2_start_time = time_ns()
    println("  Pass 2: Summing intensities into $num_bins bins...")
    
    mz_bins = range(global_min_mz, stop=global_max_mz, length=num_bins)
    bin_step = step(mz_bins)
    inv_bin_step = 1.0 / bin_step  # Precompute reciprocal
    min_mz = global_min_mz

    # Initialize thread-local intensity sums
    thread_local_intensity_sums = [zeros(Float64, num_bins) for _ in 1:Threads.nthreads()]
    thread_local_spectra_processed = [0 for _ in 1:Threads.nthreads()]

    _iterate_spectra_fast(msi_data, masked_indices === nothing ? nothing : collect(masked_indices)) do idx, mz, intensity
        thread_id = Threads.threadid()
        thread_local_spectra_processed[thread_id] += 1
        if isempty(mz)
            return
        end
        
        @inbounds @simd for i in eachindex(mz)
            # Calculate raw bin index
            raw_index = (mz[i] - min_mz) * inv_bin_step + 1.0
            bin_index = trunc(Int, raw_index)
            
            # Branchless version using clamp
            final_index = clamp(bin_index, 1, num_bins)
            thread_local_intensity_sums[thread_id][final_index] += intensity[i]
        end
    end

    # Combine thread-local results
    intensity_sum = zeros(Float64, num_bins)
    num_spectra_processed = 0
    for i in 1:Threads.nthreads()
        intensity_sum .+= thread_local_intensity_sums[i]
        num_spectra_processed += thread_local_spectra_processed[i]
    end
    pass2_duration = (time_ns() - pass2_start_time) / 1e9

    total_duration = (time_ns() - total_start_time) / 1e9
    println("\n--- mzML Profiling (Post-Optimization) ---")
    @printf "  Pass 2 (I/O + Binning):   %.2f seconds\n" pass2_duration
    @printf "  Total Function Time:      %.2f seconds\n" total_duration
    println("----------------------------------------\n")

    println("Total spectrum calculation complete for mzML.")
    return (collect(mz_bins), intensity_sum, num_spectra_processed)
end

"""
    get_total_spectrum(msi_data::MSIData; num_bins::Int=2000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the sum of all spectra in the dataset by binning.
This function dispatches to a specialized implementation based on the file type
(.imzML or .mzML) for optimal performance.

Returns a tuple containing two vectors: the binned m/z axis and the summed intensities.
"""
function get_total_spectrum(msi_data::MSIData; num_bins::Int=2000, mask_path::Union{String, Nothing}=nothing)::Tuple{Vector{Float64}, Vector{Float64}, Int}
    local masked_indices::Union{Set{Int}, Nothing} = nothing
    if mask_path !== nothing
        mask_matrix = load_and_prepare_mask(mask_path, msi_data.image_dims)
        masked_indices = get_masked_spectrum_indices(msi_data, mask_matrix)
        println("Calculating total spectrum for masked region (mask from: $(mask_path))")
    end

    if msi_data.source isa ImzMLSource
        return get_total_spectrum_imzml(msi_data, num_bins=num_bins, masked_indices=masked_indices)
    else # MzMLSource
        return get_total_spectrum_mzml(msi_data, num_bins=num_bins, masked_indices=masked_indices)
    end
end

"""
    get_average_spectrum(msi_data::MSIData; num_bins::Int=2000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the average of all spectra in the dataset by binning.
This is effectively the total spectrum divided by the number of spectra.

Returns a tuple containing two vectors: the binned m/z axis and the averaged intensities.
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
    read_compressed_array(io::IO, asset::SpectrumAsset, format::Type)

Reads a single data array (m/z or intensity) from an `.ibd` file stream,
handling both compressed and uncompressed data.

This is an internal function designed for high-performance iteration. It assumes
the file stream `io` is already positioned at the correct offset.

- If `asset.is_compressed` is true, it reads `asset.encoded_length` bytes of
  compressed data, inflates them using zlib, and reinterprets the result as a
  vector of the given `format`.
- If false, it reads `asset.encoded_length` *elements* of uncompressed data
  directly into a vector.

# Arguments
- `io`: The IO stream of the `.ibd` file.
- `asset`: The `SpectrumAsset` for the array.
- `format`: The data type of the elements in the array.

# Returns
- A `Vector` containing the data.

# Throws
- An error if zlib decompression fails, which can indicate corrupt data or
  an incorrect offset in the `.imzML` metadata.
"""
function read_compressed_array(data::MSIData, io::IO, asset::SpectrumAsset, format::Type)
    # Add validation before seeking
    if asset.offset < 0 || asset.offset >= filesize(io)
        throw(FileFormatError("Invalid asset offset: $(asset.offset) for file size $(filesize(io))"))
    end
    seek(io, asset.offset)
    
    if asset.is_compressed
        # Get buffer for compressed bytes
        compressed_bytes_buffer = get_buffer!(data.buffer_pool, asset.encoded_length)
        readbytes!(io, compressed_bytes_buffer, asset.encoded_length)
        
        println("DEBUG: Decompressing data - offset=$(asset.offset), compressed_bytes=$(length(compressed_bytes_buffer))")
        
        local decompressed_bytes_buffer
        try
            # Estimate decompressed size (can be larger than compressed)
            # A common heuristic is 4x compressed size, but zlib can be more efficient
            # For now, let Libz.inflate handle allocation, then copy to pooled buffer
            # This is a temporary allocation, will be optimized later if needed
            temp_decompressed = Libz.inflate(compressed_bytes_buffer)
            
            decompressed_bytes_buffer = get_buffer!(data.buffer_pool, length(temp_decompressed))
            copyto!(decompressed_bytes_buffer, temp_decompressed)

            println("DEBUG: Decompression successful - decompressed_bytes=$(length(decompressed_bytes_buffer))")
        catch e
            @error "ZLIB DECOMPRESSION FAILED. This is likely due to an incorrect offset or corrupt data in the .ibd file."
            @error "Asset offset: $(asset.offset), Encoded length: $(asset.encoded_length)"
            # Print first 16 bytes to stderr for diagnosis
            bytes_to_print = min(16, length(compressed_bytes_buffer))
            @error "First $bytes_to_print bytes of the data chunk we tried to decompress:"
            println(stderr, view(compressed_bytes_buffer, 1:bytes_to_print))
            rethrow(e)
        finally
            release_buffer!(data.buffer_pool, compressed_bytes_buffer)
        end
        
        # Use an IOBuffer to safely read the data
        bytes_io = IOBuffer(decompressed_bytes_buffer)
        n_elements = bytes_io.size ÷ sizeof(format)
        array = Array{format}(undef, n_elements)
        read!(bytes_io, array)
        
        release_buffer!(data.buffer_pool, decompressed_bytes_buffer)
        return array
    else
        # Read uncompressed data directly
        # For uncompressed imzML, encoded_length is the number of elements
        array = Vector{format}(undef, asset.encoded_length)
        read!(io, array)
        return array
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
"""
function _iterate_uncompressed_fast(f::Function, data::MSIData, source::ImzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    # Optimized path for uncompressed data using buffer reuse

    # Determine which indices to iterate over
    spectrum_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : collect(indices_to_iterate) # Ensure it's a concrete vector

    # Determine max_points for pre-allocation (common for both paths)
    max_points = isempty(data.spectra_metadata) ? 0 : maximum(meta -> meta.mz_asset.encoded_length, data.spectra_metadata)
    
    # Pre-allocate buffers for the output (mutable arrays for in-place ltoh)
    mz_buffer = Vector{source.mz_format}(undef, max_points)
    int_buffer = Vector{source.intensity_format}(undef, max_points)

    for i in spectrum_indices
        meta = data.spectra_metadata[i]
        nPoints = meta.mz_asset.encoded_length

        if nPoints == 0
            f(i, view(mz_buffer, 0:-1), view(int_buffer, 0:-1))
            continue
        end

        mz_view = view(mz_buffer, 1:nPoints)
        int_view = view(int_buffer, 1:nPoints)

        # Use appropriate read mechanism
        if source.ibd_handle isa ThreadSafeFileHandle
            # Ensure correct order of reading if offsets are close
            if meta.mz_asset.offset < meta.int_asset.offset
                read_at!(source.ibd_handle, mz_view, meta.mz_asset.offset)
                read_at!(source.ibd_handle, int_view, meta.int_asset.offset)
            else
                read_at!(source.ibd_handle, int_view, meta.int_asset.offset)
                read_at!(source.ibd_handle, mz_view, meta.mz_asset.offset)
            end
        else # Plain IOStream
            # Ensure correct order of reading if offsets are close
            if meta.mz_asset.offset < meta.int_asset.offset
                seek(source.ibd_handle, meta.mz_asset.offset)
                read!(source.ibd_handle, mz_view)
                seek(source.ibd_handle, meta.int_asset.offset) # FIX: Added missing seek
                read!(source.ibd_handle, int_view)
            else
                seek(source.ibd_handle, meta.int_asset.offset)
                read!(source.ibd_handle, int_view)
                seek(source.ibd_handle, meta.mz_asset.offset) # FIX: Added missing seek
                read!(source.ibd_handle, mz_view)
            end
        end

        mz_view .= ltoh.(mz_view)
        int_view .= ltoh.(int_view)
        f(i, mz_view, int_view)
    end
end

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
"""
function _iterate_compressed_fast(f::Function, data::MSIData, source::ImzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    # Path for datasets containing at least one compressed spectrum.
    # This path reads and decompresses each spectrum individually.
    
    # Determine which indices to iterate over
    spectrum_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate

    for i in spectrum_indices
        meta = data.spectra_metadata[i]
        
        if meta.mz_asset.encoded_length == 0 && meta.int_asset.encoded_length == 0
            f(i, source.mz_format[], source.intensity_format[])
            continue
        end

        # Read and decompress each array
        mz_array = read_compressed_array(data, source.ibd_handle, meta.mz_asset, source.mz_format)
        intensity_array = read_compressed_array(data, source.ibd_handle, meta.int_asset, source.intensity_format)
        
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
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    if isempty(data.spectra_metadata)
        return
    end
    
    # Check if ANY spectra are compressed and dispatch to the appropriate implementation
    any_compressed = any(meta -> meta.mz_asset.is_compressed || meta.int_asset.is_compressed, 
                        data.spectra_metadata)
    
    if any_compressed
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
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::MzMLSource, indices_to_iterate::Union{AbstractVector{Int}, Nothing})
    # This implementation is for mzML. To improve disk I/O, we can reorder the read
    # operations to be as sequential as possible based on their offset in the file.
    
    # Determine which indices to iterate over
    spectrum_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate

    # Create a vector of (index, offset) tuples to be sorted
    indices_with_offsets = [(i, data.spectra_metadata[i].mz_asset.offset) for i in spectrum_indices]
    
    # Sort by offset to make disk access more sequential
    sort!(indices_with_offsets, by = x -> x[2])

    for (i, _) in indices_with_offsets
        meta = data.spectra_metadata[i]
        
        mz = read_binary_vector(data, source.file_handle, meta.mz_asset)
        intensity = read_binary_vector(data, source.file_handle, meta.int_asset)
        
        f(i, mz, intensity)
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
number of available threads (`Threads.nthreads()`).

# Arguments
- `f`: A function to call for each spectrum, with the signature `f(index, mz, intensity)`.
- `data`: The `MSIData` object.
- `indices_to_iterate`: An optional `AbstractVector{Int}` specifying a subset of
  spectrum indices to iterate over. If `nothing`, it iterates over all spectra.
"""
function _iterate_spectra_fast(f::Function, data::MSIData, indices_to_iterate::Union{AbstractVector{Int}, Nothing}=nothing)
    # Default to all indices if not specified
    all_indices = (indices_to_iterate === nothing) ? (1:length(data.spectra_metadata)) : indices_to_iterate
    
    if Threads.nthreads() > 1 && length(all_indices) > 100 # Heuristic: only parallelize for enough work
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
    # Split indices into chunks for each thread
    n_chunks = Threads.nthreads()
    chunk_size = ceil(Int, length(indices) / n_chunks)
    chunks = collect(Iterators.partition(indices, chunk_size))
    
    Threads.@threads for chunk in chunks
        # Each thread gets its own file handle based on source type
        if data.source isa ImzMLSource
            local_handle = open(data.source.ibd_handle.path, "r")
            local_source = ImzMLSource(local_handle, data.source.mz_format, data.source.intensity_format)
            
            try
                # Use the appropriate implementation with thread-local source
                _iterate_spectra_fast_impl(f, data, local_source, chunk)
            finally
                close(local_handle)
            end
        elseif data.source isa MzMLSource
            local_handle = open(data.source.file_handle.path, "r")
            local_source = MzMLSource(local_handle, data.source.mz_format, data.source.intensity_format)
            
            try
                _iterate_spectra_fast_impl(f, data, local_source, chunk)
            finally
                close(local_handle)
            end
        end
    end
end

"""
    _iterate_spectra_fast_parallel(f::Function, data::MSIData)

Parallel version that processes all spectra in the dataset.
Convenience wrapper for the indexed version.
"""
function _iterate_spectra_fast_parallel(f::Function, data::MSIData)
    all_indices = collect(1:length(data.spectra_metadata))
    _iterate_spectra_fast_parallel(f, data, all_indices)
end

"""
    _iterate_spectra_fast_threadsafe(f::Function, data::MSIData)

Thread-safe version that processes all spectra in the dataset.
Convenience wrapper for the indexed version.
"""
function _iterate_spectra_fast_threadsafe(f::Function, data::MSIData)
    all_indices = collect(1:length(data.spectra_metadata))
    _iterate_spectra_fast_threadsafe(f, data, all_indices)
end

"""
    get_bloom_filter(data::MSIData, index::Int)::BloomFilter{Int}

Retrieves the Bloom filter for a specific spectrum by its index.
If the Bloom filter has not yet been created (due to lazy initialization),
its will be created on-demand and cached.

# Arguments
- `data`: The `MSIData` object.
- `index`: The linear index of the spectrum.

# Returns
- A `BloomFilter{Int}` for the specified spectrum.
"""
function get_bloom_filter(data::MSIData, index::Int)::BloomFilter{Int}
    lock(data.cache_lock) do
        if data.bloom_filters === nothing
            # Initialize Bloom filters on-demand
            data.bloom_filters = Vector{Union{BloomFilter{Int}, Nothing}}(undef, length(data.spectra_metadata))
            fill!(data.bloom_filters, nothing) # Initialize with nothing
        end
        
        if data.bloom_filters[index] === nothing
            # Create Bloom filter lazily
            mz, intensity = GetSpectrum(data, index)
            data.bloom_filters[index] = create_bloom_filter_for_spectrum(mz)
        end
        
        return data.bloom_filters[index]
    end
end
