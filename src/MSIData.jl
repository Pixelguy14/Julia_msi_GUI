# src/MSIData.jl

"""
This file defines the unified MSIData object and the associated data access layer,
including caching and iteration logic, for handling large mzML and imzML datasets
efficiently.
"""

using Base64, Libz, Serialization, Printf # For reading binary data

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
    ibd_handle::IO
    mz_format::Type
    intensity_format::Type
end

"""
    MzMLSource <: MSDataSource

A data source for `.mzML` files, holding a handle to the `.mzML` file itself
(which contains the binary data encoded in Base64) and the expected data formats.
"""
struct MzMLSource <: MSDataSource
    file_handle::IO
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
- `encoded_length`: The length of the data. For mzML, this is the Base64 encoded length.
  For imzML, this is the number of elements in the array.
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
    SpectrumMetadata

Contains all metadata for a single spectrum, common to both imzML and mzML formats.

# Fields
- `x`, `y`: The spatial coordinates of the spectrum (for imzML only).
- `id`: The unique identifier string for the spectrum (for mzML only).
- `mz_asset`: A `SpectrumAsset` for the m/z array.
- `int_asset`: A `SpectrumAsset` for the intensity array.
"""
struct SpectrumMetadata
    # For imzML
    x::Int32
    y::Int32
    
    # For mzML
    id::String

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
- `image_dims`: A tuple `(width, height)` of the spatial dimensions (for imzML).
- `coordinate_map`: A matrix mapping `(x, y)` coordinates to a linear spectrum index (for imzML).
- `cache`: A dictionary holding cached spectra.
- `cache_order`: A vector tracking the usage order for the LRU cache.
- `cache_size`: The maximum number of spectra to store in the cache.
"""
mutable struct MSIData
    source::MSDataSource
    spectra_metadata::Vector{SpectrumMetadata}
    image_dims::Tuple{Int, Int} # (width, height) for imaging data
    coordinate_map::Union{Matrix{Int}, Nothing} # Maps (x,y) to linear index for imzML

    # LRU Cache implementation
    cache::Dict{Int, Tuple{Vector, Vector}}
    cache_order::Vector{Int} # Stores indices, with most recently used at the end
    cache_size::Int # Max number of spectra in cache

    function MSIData(source, metadata, dims, coordinate_map, cache_size)
        obj = new(source, metadata, dims, coordinate_map, Dict(), [], cache_size)
        
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


# --- Internal function for reading binary data --- #

"""
    read_binary_vector(io::IO, asset::SpectrumAsset)

Reads and decodes a single binary data vector (like m/z or intensity array)
from a `.mzML` file. The data is expected to be Base64-encoded and may be
compressed.

This internal function handles:
1. Reading the raw Base64 string.
2. Decoding from Base64.
3. Decompressing the data if `asset.is_compressed` is true.
4. Converting the byte order from network (big-endian) to host order.

# Arguments
- `io`: The IO stream of the `.mzML` file.
- `asset`: The `SpectrumAsset` containing metadata for the array.

# Returns
- A `Vector` of the appropriate type containing the decoded data.
"""
function read_binary_vector(io::IO, asset::SpectrumAsset)
    seek(io, asset.offset)
    raw_b64 = read(io, asset.encoded_length)
    decoded_bytes = base64decode(strip(String(raw_b64)))
    bytes_io = IOBuffer(asset.is_compressed ? Libz.inflate(decoded_bytes) : decoded_bytes)
    
    n_elements = Int(bytes_io.size / sizeof(asset.format))
    out_array = Array{asset.format}(undef, n_elements)
    read!(bytes_io, out_array)
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
function read_spectrum_from_disk(source::ImzMLSource, meta::SpectrumMetadata)
    # For imzML, the binary data is raw, not base64 encoded. 
    # The `encoded_length` field in this case holds the number of points.
    mz = Array{source.mz_format}(undef, meta.mz_asset.encoded_length)
    intensity = Array{source.intensity_format}(undef, meta.int_asset.encoded_length)

    seek(source.ibd_handle, meta.mz_asset.offset)
    read!(source.ibd_handle, mz)

    seek(source.ibd_handle, meta.int_asset.offset)
    read!(source.ibd_handle, intensity)

    # imzML data is little-endian. Convert to host byte order.
    mz .= ltoh.(mz)
    intensity .= ltoh.(intensity)

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
function read_spectrum_from_disk(source::MzMLSource, meta::SpectrumMetadata)
    # For mzML, data is Base64 encoded within the XML
    mz = read_binary_vector(source.file_handle, meta.mz_asset)
    intensity = read_binary_vector(source.file_handle, meta.int_asset)
    return mz, intensity
end

# --- Public API --- #

"""
    GetSpectrum(data::MSIData, index::Int)

Retrieves a single spectrum by its index, utilizing a cache for performance.

This function is the core of the "Indexed" and "Cache" access patterns.
"""
function GetSpectrum(data::MSIData, index::Int)
    if index < 1 || index > length(data.spectra_metadata)
        error("Spectrum index $index out of bounds.")
    end

    # Phase 1: Check the cache
    if haskey(data.cache, index)
        # Cache Hit: Move item to the end of the LRU list and return from cache
        filter!(x -> x != index, data.cache_order)
        push!(data.cache_order, index)
        return data.cache[index]
    end

    # Phase 2: Cache Miss - Read from disk
    meta = data.spectra_metadata[index]
    spectrum = read_spectrum_from_disk(data.source, meta)

    # Phase 3: Update cache
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

"""
    GetSpectrum(data::MSIData, x::Int, y::Int)

Retrieves a single spectrum by its (x, y) coordinates for imaging data.
Utilizes a coordinate map for efficient lookup and then the cache.
"""
function GetSpectrum(data::MSIData, x::Int, y::Int)
    if data.coordinate_map === nothing
        error("Coordinate map not available. This method is only for imaging data loaded from .imzML files.")
    end
    width, height = data.image_dims
    if x < 1 || x > width || y < 1 || y > height
        error("Coordinates ($x, $y) out of bounds for image dimensions ($width, $height).")
    end

    index = data.coordinate_map[x, y]
    if index == 0
        error("No spectrum found at coordinates ($x, $y).")
    end

    return GetSpectrum(data, index) # Call the existing method
end

using Serialization

function get_total_spectrum_imzml(msi_data::MSIData; num_bins::Int=2000)
    println("Calculating total spectrum for imzML (2-pass method)...")
    total_start_time = time_ns()

    # 1. First Pass: Find the global m/z range
    pass1_start_time = time_ns()
    println("  Pass 1: Finding global m/z range...")
    global_min_mz = Inf
    global_max_mz = -Inf
    _iterate_spectra_fast(msi_data) do idx, mz, _
        if !isempty(mz)
            local_min, local_max = extrema(mz)
            global_min_mz = min(global_min_mz, local_min)
            global_max_mz = max(global_max_mz, local_max)
        end
    end
    pass1_duration = (time_ns() - pass1_start_time) / 1e9

    if !isfinite(global_min_mz)
        @warn "Could not determine a valid m/z range for imzML. All spectra might be empty."
        return (Float64[], Float64[])
    end
    println("  Global m/z range found: [$(global_min_mz), $(global_max_mz)]")

    # 2. Define Bins and precompute constants
    mz_bins = range(global_min_mz, stop=global_max_mz, length=num_bins)
    intensity_sum = zeros(Float64, num_bins)
    bin_step = step(mz_bins)
    inv_bin_step = 1.0 / bin_step  # Precompute reciprocal to avoid division
    min_mz = global_min_mz

    # 3. Second Pass: Optimized binning
    pass2_start_time = time_ns()
    println("  Pass 2: Summing intensities into $num_bins bins...")
    _iterate_spectra_fast(msi_data) do idx, mz, intensity
        if isempty(mz)
            return
        end
        
        # Use views to avoid bounds checks on the input arrays
        mz_view = mz
        intensity_view = intensity
        
        # Manual binning without clamp/round for SIMD
        @inbounds for i in eachindex(mz_view)
            # Calculate raw bin index (much faster than round+clamp)
            raw_index = (mz_view[i] - min_mz) * inv_bin_step + 1.0
            bin_index = trunc(Int, raw_index)
            
            # Manual bounds checking (faster than clamp)
            if 1 <= bin_index <= num_bins
                intensity_sum[bin_index] += intensity_view[i]
            elseif bin_index < 1
                intensity_sum[1] += intensity_view[i]
            else # bin_index > num_bins
                intensity_sum[num_bins] += intensity_view[i]
            end
        end
    end
    pass2_duration = (time_ns() - pass2_start_time) / 1e9

    total_duration = (time_ns() - total_start_time) / 1e9
    println("\n--- imzML Profiling ---")
    @printf "  Pass 1 (I/O only):        %.2f seconds\n" pass1_duration
    @printf "  Pass 2 (I/O + Binning):   %.2f seconds\n" pass2_duration
    @printf "  Est. Binning Overhead:    %.2f seconds\n" (pass2_duration - pass1_duration)
    @printf "  Total Function Time:      %.2f seconds\n" total_duration
    println("-------------------------\n")

    println("Total spectrum calculation complete for imzML.")
    return (collect(mz_bins), intensity_sum)
end

function get_total_spectrum_mzml(msi_data::MSIData; num_bins::Int=2000)
    println("Calculating total spectrum for mzML (optimized single-pass method)...")
    total_start_time = time_ns()
    num_spectra = length(msi_data.spectra_metadata)
    if num_spectra == 0
        return (Float64[], Float64[])
    end

    temp_path, temp_io = mktemp()
    try
        # --- Pass 1: Read from source, find m/z range, and write decoded spectra to temp file ---
        pass1_start_time = time_ns()
        println("  Pass 1: Caching decoded spectra and finding global m/z range...")
        global_min_mz = Inf
        global_max_mz = -Inf

        _iterate_spectra_fast(msi_data) do idx, mz, intensity
            if !isempty(mz)
                local_min, local_max = extrema(mz)
                global_min_mz = min(global_min_mz, local_min)
                global_max_mz = max(global_max_mz, local_max)
            end
            Serialization.serialize(temp_io, (mz, intensity))
        end

        flush(temp_io)
        pass1_duration = (time_ns() - pass1_start_time) / 1e9

        if !isfinite(global_min_mz)
            @warn "Could not determine a valid m/z range for mzML. All spectra might be empty."
            return (Float64[], Float64[])
        end
        println("  Global m/z range found: [$(global_min_mz), $(global_max_mz)]")

        # --- Pass 2: Read from fast temp file and bin intensities ---
        pass2_start_time = time_ns()
        println("  Pass 2: Reading from cache and summing intensities into $num_bins bins...")
        seekstart(temp_io)
        mz_bins = range(global_min_mz, stop=global_max_mz, length=num_bins)
        intensity_sum = zeros(Float64, num_bins)
        bin_step = step(mz_bins)

        while !eof(temp_io)
            mz, intensity = Serialization.deserialize(temp_io)::Tuple{AbstractVector, AbstractVector}
            if isempty(mz)
                continue
            end
            for i in eachindex(mz)
                bin_index = clamp(round(Int, (mz[i] - global_min_mz) / bin_step) + 1, 1, num_bins)
                intensity_sum[bin_index] += intensity[i]
            end
        end
        pass2_duration = (time_ns() - pass2_start_time) / 1e9

        total_duration = (time_ns() - total_start_time) / 1e9
        println("\n--- mzML Profiling ---")
        @printf "  Pass 1 (Read+Decode+Cache): %.2f seconds\n" pass1_duration
        @printf "  Pass 2 (Read Cache+Bin):    %.2f seconds\n" pass2_duration
        @printf "  Total Function Time:        %.2f seconds\n" total_duration
        println("----------------------\n")

        println("Total spectrum calculation complete for mzML.")
        return (collect(mz_bins), intensity_sum)

    finally
        close(temp_io)
        rm(temp_path, force=true)
        println("  Temporary cache file removed.")
    end
end

"""
    get_total_spectrum(msi_data::MSIData; num_bins::Int=2000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the sum of all spectra in the dataset by binning.
This function dispatches to a specialized implementation based on the file type
(.imzML or .mzML) for optimal performance.

Returns a tuple containing two vectors: the binned m/z axis and the summed intensities.
"""
function get_total_spectrum(msi_data::MSIData; num_bins::Int=2000)
    if msi_data.source isa ImzMLSource
        return get_total_spectrum_imzml(msi_data, num_bins=num_bins)
    else # MzMLSource
        return get_total_spectrum_mzml(msi_data, num_bins=num_bins)
    end
end

"""
    get_average_spectrum(msi_data::MSIData; num_bins::Int=20000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the average of all spectra in the dataset by binning.
This is effectively the total ion chromatogram (TIC) divided by the number of spectra.

Returns a tuple containing two vectors: the binned m/z axis and the averaged intensities.
"""
function get_average_spectrum(msi_data::MSIData; num_bins::Int=2000)
    # This function uses the exact same logic as get_total_spectrum...
    mz_bins, intensity_sum = get_total_spectrum(msi_data, num_bins=num_bins)

    if isempty(intensity_sum)
        return (mz_bins, intensity_sum)
    end

    # ...with one final step: dividing by the number of spectra to get the average.
    println("Averaging spectrum...")
    num_spectra = length(msi_data.spectra_metadata)
    average_intensity = intensity_sum ./ num_spectra
    
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

Returns an iterator that yields each spectrum, processing the file sequentially
with minimal memory overhead. This iterator supports caching via `GetSpectrum`.

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
    _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource)

Internal implementation of the fast iterator for `.imzML` files. It reads
data directly from the `.ibd` file stream and reuses pre-allocated buffers
to minimize memory allocations and overhead, making it ideal for bulk processing tasks.
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource)
    # Optimized implementation: allocates arrays per spectrum and uses an efficient I/O pattern.
    for i in 1:length(data.spectra_metadata)
        meta = data.spectra_metadata[i]
        nPoints = meta.mz_asset.encoded_length

        if nPoints == 0
            f(i, source.mz_format[], source.intensity_format[])
            continue
        end

        # Allocate fresh arrays for each spectrum.
        mz = Vector{source.mz_format}(undef, nPoints)
        intensity = Vector{source.intensity_format}(undef, nPoints)

        # --- I/O OPTIMIZATION: Seek only once per spectrum ---
        # The mz and intensity data are contiguous, so after reading the first,
        # we can immediately read the second without a costly second seek.
        if meta.mz_asset.offset < meta.int_asset.offset
            # m/z is first, seek to it.
            seek(source.ibd_handle, meta.mz_asset.offset)
            # Read m/z, then immediately read intensity.
            read!(source.ibd_handle, mz)
            read!(source.ibd_handle, intensity)
        else
            # intensity is first, seek to it.
            seek(source.ibd_handle, meta.int_asset.offset)
            # Read intensity, then immediately read m/z.
            read!(source.ibd_handle, intensity)
            read!(source.ibd_handle, mz)
        end

        # Convert byte order.
        mz .= ltoh.(mz)
        intensity .= ltoh.(intensity)

        f(i, mz, intensity)
    end
end

"""
    _iterate_spectra_fast_impl(f::Function, data::MSIData, source::MzMLSource)

Internal implementation of the fast iterator for `.mzML` files. It iterates
through each spectrum, decodes the Base64 data on the fly, and calls the
provided function. It bypasses the `GetSpectrum` cache to avoid storing
all decoded spectra in memory.
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::MzMLSource)
    # This implementation is for mzML. It iterates through each spectrum,
    # decodes it, and then calls the function `f`. It's less performant than
    # the imzML version because we cannot pre-allocate buffers of a known size,
    # but it is correct and still faster than using GetSpectrum due to bypassing the cache.
    for i in 1:length(data.spectra_metadata)
        meta = data.spectra_metadata[i]
        
        mz = read_binary_vector(source.file_handle, meta.mz_asset)
        intensity = read_binary_vector(source.file_handle, meta.int_asset)
        
        f(i, mz, intensity)
    end
end

"""
    _iterate_spectra_fast(f::Function, data::MSIData)

An internal, high-performance iterator for bulk processing that bypasses the cache.
It dispatches to a specialized implementation based on the data source type
(`ImzMLSource` or `MzMLSource`).

- `f`: A function to call for each spectrum, with signature `f(index, mz, intensity)`.
- `data`: The `MSIData` object.
"""
function _iterate_spectra_fast(f::Function, data::MSIData)
    # Dispatch to the correct implementation based on the source type
    _iterate_spectra_fast_impl(f, data, data.source)
end
