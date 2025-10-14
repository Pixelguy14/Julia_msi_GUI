# src/MSIData.jl

"""
This file defines the unified MSIData object and the associated data access layer,
including caching and iteration logic, for handling large mzML and imzML datasets
efficiently.
"""

using Base64, Libz, Serialization, Printf, DataFrames, Base.Threads # For reading binary data

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
- `encoded_length`: The length of the data. For uncompressed imzML, this is the number of 
  elements in the array. For compressed imzML, it is the number of bytes of the compressed 
  data. For mzML, this is the length of the Base64 encoded string.
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

@enum SpectrumMode CENTROID=1 PROFILE=2 UNKNOWN=3

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
    image_dims::Tuple{Int, Int} # (width, height) for imaging data
    coordinate_map::Union{Matrix{Int}, Nothing} # Maps (x,y) to linear index for imzML

    # LRU Cache for GetSpectrum
    cache::Dict{Int, Tuple{Vector, Vector}}
    cache_order::Vector{Int} # Stores indices, with most recently used at the end
    cache_size::Int # Max number of spectra in cache
    cache_lock::ReentrantLock # To make cache access thread-safe

    # Pre-computed analytics/metadata
    global_min_mz::Union{Float64, Nothing}
    global_max_mz::Union{Float64, Nothing}
    spectrum_stats_df::Union{DataFrame, Nothing}

    function MSIData(source, metadata, dims, coordinate_map, cache_size)
        obj = new(source, metadata, dims, coordinate_map, 
                  Dict(), [], cache_size, ReentrantLock(),
                  nothing, nothing, nothing) # Initialize new fields to nothing
        
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
        error("Spectrum index $index out of bounds.")
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
    spectrum = read_spectrum_from_disk(data.source, meta)

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
    # Idempotency check
    if msi_data.spectrum_stats_df !== nothing && hasproperty(msi_data.spectrum_stats_df, :MinMZ)
        println("Analytics have already been pre-computed.")
        return
    end
    
    println("Pre-computing analytics (single pass)...")
    start_time = time_ns()

    num_spectra = length(msi_data.spectra_metadata)
    
    # DEBUG: Check the first spectrum's metadata
    if num_spectra > 0
        first_meta = msi_data.spectra_metadata[1]
        println("DEBUG First spectrum metadata:")
        println("  mz_asset: format=$(first_meta.mz_asset.format), compressed=$(first_meta.mz_asset.is_compressed)")
        println("  mz_asset: offset=$(first_meta.mz_asset.offset), encoded_length=$(first_meta.mz_asset.encoded_length)")
        println("  int_asset: format=$(first_meta.int_asset.format), compressed=$(first_meta.int_asset.is_compressed)")
        println("  int_asset: offset=$(first_meta.int_asset.offset), encoded_length=$(first_meta.int_asset.encoded_length)")
        println("  mode: $(first_meta.mode)")
    end

    # Initialize variables for global stats
    g_min_mz = Inf
    g_max_mz = -Inf

    # Initialize vectors for per-spectrum stats
    tics = Vector{Float64}(undef, num_spectra)
    bpis = Vector{Float64}(undef, num_spectra)
    bp_mzs = Vector{Float64}(undef, num_spectra)
    num_points = Vector{Int}(undef, num_spectra)
    min_mzs = Vector{Float64}(undef, num_spectra)
    max_mzs = Vector{Float64}(undef, num_spectra)
    modes = Vector{SpectrumMode}(undef, num_spectra)
    is_compressed = Vector{Bool}(undef, num_spectra)

    # DEBUG: Add counter to see how many spectra have data
    spectra_with_data = 0
    empty_spectra = 0

    _iterate_spectra_fast(msi_data) do idx, mz, intensity
        # Store metadata
        modes[idx] = msi_data.spectra_metadata[idx].mode
        is_compressed[idx] = msi_data.spectra_metadata[idx].mz_asset.is_compressed || 
                            msi_data.spectra_metadata[idx].int_asset.is_compressed
        
        if isempty(mz)
            empty_spectra += 1
            tics[idx] = 0.0
            bpis[idx] = 0.0
            bp_mzs[idx] = 0.0
            num_points[idx] = 0
            min_mzs[idx] = Inf
            max_mzs[idx] = -Inf
            return
        else
            spectra_with_data += 1
        end

        # Update global m/z range
        local_min, local_max = extrema(mz)
        g_min_mz = min(g_min_mz, local_min)
        g_max_mz = max(g_max_mz, local_max)
        min_mzs[idx] = local_min
        max_mzs[idx] = local_max

        # Calculate per-spectrum stats
        tics[idx] = sum(intensity)
        max_int, max_idx = findmax(intensity)
        bpis[idx] = max_int
        bp_mzs[idx] = mz[max_idx]
        num_points[idx] = length(mz)
    end

    # Add mode statistics
    centroid_count = count(==(CENTROID), modes)
    profile_count = count(==(PROFILE), modes)
    unknown_count = count(==(UNKNOWN), modes)
    
    println("DEBUG Mode Statistics:")
    println("  Centroid spectra: $centroid_count")
    println("  Profile spectra: $profile_count")
    println("  Unknown mode: $unknown_count")

    # DEBUG: Print summary
    println("DEBUG Analytics Summary:")
    println("  Total spectra: $num_spectra")
    println("  Spectra with data: $spectra_with_data")
    println("  Empty spectra: $empty_spectra")
    println("  Global m/z range: [$g_min_mz, $g_max_mz]")
    println("  Centroid spectra: $(count(==(CENTROID), modes))")
    println("  Profile spectra: $(count(==(PROFILE), modes))")
    println("  Compressed spectra: $(sum(is_compressed))")
    println("  Average points per spectrum: $(mean(num_points))")

    # Populate the MSIData object
    msi_data.global_min_mz = g_min_mz
    msi_data.global_max_mz = g_max_mz
    msi_data.spectrum_stats_df = DataFrame(
        SpectrumID = 1:num_spectra,
        TIC = tics,
        BPI = bpis,
        BasePeakMZ = bp_mzs,
        NumPoints = num_points,
        MinMZ = min_mzs,
        MaxMZ = max_mzs,
        Mode = modes,
        IsCompressed = is_compressed
    )
    
    duration = (time_ns() - start_time) / 1e9
    @printf "Analytics pre-computation complete in %.2f seconds.\n" duration
    
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
function get_total_spectrum_imzml(msi_data::MSIData; num_bins::Int=2000)
    println("Calculating total spectrum for imzML (2-pass method)...")
    total_start_time = time_ns()

    local global_min_mz, global_max_mz

    if msi_data.global_min_mz !== nothing
        println("  Using pre-computed m/z range.")
        global_min_mz = msi_data.global_min_mz
        global_max_mz = msi_data.global_max_mz
    else
        # 1. First Pass: Find the global m/z range by reading the fast .ibd file
        pass1_start_time = time_ns()
        println("  Pass 1: Finding global m/z range...")
        g_min_mz = Inf
        g_max_mz = -Inf
        _iterate_spectra_fast(msi_data) do idx, mz, _
            if !isempty(mz)
                local_min, local_max = extrema(mz)
                g_min_mz = min(g_min_mz, local_min)
                g_max_mz = max(g_max_mz, local_max)
            end
        end
        pass1_duration = (time_ns() - pass1_start_time) / 1e9
        if !isfinite(g_min_mz)
            @warn "Could not determine a valid m/z range for imzML. All spectra might be empty."
            return (Float64[], Float64[])
        end
        
        # Use and cache the result
        global_min_mz = g_min_mz
        global_max_mz = g_max_mz
        msi_data.global_min_mz = global_min_mz
        msi_data.global_max_mz = global_max_mz
        println("  Global m/z range found and cached: [$(global_min_mz), $(global_max_mz)]")
        @printf "  (Pass 1 took %.2f seconds)\n" pass1_duration
    end
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
    println("\n--- imzML Profiling (Post-Optimization) ---")
    @printf "  Pass 2 (I/O + Binning):   %.2f seconds\n" pass2_duration
    @printf "  Total Function Time:      %.2f seconds\n" total_duration
    println("----------------------------------------\n")

    println("Total spectrum calculation complete for imzML.")
    return (collect(mz_bins), intensity_sum)
end

"""
    get_total_spectrum_mzml(msi_data::MSIData; num_bins::Int=2000)

Internal function to calculate the total spectrum for an `.mzML` file.

It uses a two-pass approach analogous to the `imzML` implementation:
1. The first pass finds the global m/z range by iterating through all spectra.
2. The second pass sums intensities into a pre-defined number of bins.

This function is called by `get_total_spectrum`.
"""
function get_total_spectrum_mzml(msi_data::MSIData; num_bins::Int=2000)
    println("Calculating total spectrum for mzML (2-pass method)...")
    total_start_time = time_ns()

    local global_min_mz, global_max_mz

    if msi_data.global_min_mz !== nothing
        println("  Using pre-computed m/z range.")
        global_min_mz = msi_data.global_min_mz
        global_max_mz = msi_data.global_max_mz
    else
        # --- Pass 1: Find m/z range and cache it ---
        pass1_start_time = time_ns()
        println("  Pass 1: Finding global m/z range...")
        g_min_mz = Inf
        g_max_mz = -Inf
        _iterate_spectra_fast(msi_data) do idx, mz, intensity
            if !isempty(mz)
                local_min, local_max = extrema(mz)
                g_min_mz = min(g_min_mz, local_min)
                g_max_mz = max(g_max_mz, local_max)
            end
        end
        pass1_duration = (time_ns() - pass1_start_time) / 1e9

        if !isfinite(g_min_mz)
            @warn "Could not determine a valid m/z range for mzML. All spectra might be empty."
            return (Float64[], Float64[])
        end

        # Use and cache the result
        global_min_mz = g_min_mz
        global_max_mz = g_max_mz
        msi_data.global_min_mz = global_min_mz
        msi_data.global_max_mz = global_max_mz

        println("  Global m/z range found and cached: [$(global_min_mz), $(global_max_mz)]")
        @printf "  (Pass 1 took %.2f seconds)\n" pass1_duration
    end

    # --- Pass 2: Bin intensities ---
    pass2_start_time = time_ns()
    println("  Pass 2: Summing intensities into $num_bins bins...")
    
    mz_bins = range(global_min_mz, stop=global_max_mz, length=num_bins)
    intensity_sum = zeros(Float64, num_bins)
    bin_step = step(mz_bins)
    inv_bin_step = 1.0 / bin_step  # Precompute reciprocal
    min_mz = global_min_mz

    _iterate_spectra_fast(msi_data) do idx, mz, intensity
        if isempty(mz)
            return
        end
        
        @inbounds for i in eachindex(mz)
            # Calculate raw bin index
            raw_index = (mz[i] - min_mz) * inv_bin_step + 1.0
            bin_index = trunc(Int, raw_index)
            
            # Manual bounds checking
            if 1 <= bin_index <= num_bins
                intensity_sum[bin_index] += intensity[i]
            elseif bin_index < 1
                intensity_sum[1] += intensity[i]
            else # bin_index > num_bins
                intensity_sum[num_bins] += intensity[i]
            end
        end
    end
    pass2_duration = (time_ns() - pass2_start_time) / 1e9

    total_duration = (time_ns() - total_start_time) / 1e9
    println("\n--- mzML Profiling (Post-Optimization) ---")
    @printf "  Pass 2 (I/O + Binning):   %.2f seconds\n" pass2_duration
    @printf "  Total Function Time:      %.2f seconds\n" total_duration
    println("----------------------------------------\n")

    println("Total spectrum calculation complete for mzML.")
    return (collect(mz_bins), intensity_sum)
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
    get_average_spectrum(msi_data::MSIData; num_bins::Int=2000) -> Tuple{Vector{Float64}, Vector{Float64}}

Calculates the average of all spectra in the dataset by binning.
This is effectively the total spectrum divided by the number of spectra.

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

- If `asset.is_compressed` is true, it reads the compressed bytes, inflates
  them using zlib, and reinterprets the result as a vector of the given `format`.
- If false, it reads the uncompressed data directly into a vector.

# Arguments
- `io`: The IO stream of the `.ibd` file.
- `asset`: The `SpectrumAsset` for the array.
- `format`: The data type of the elements in the array.

# Returns
- A `Vector` containing the data.
"""
function read_compressed_array(io::IO, asset::SpectrumAsset, format::Type)
    seek(io, asset.offset)
    
    if asset.is_compressed
        # Read compressed bytes - use encoded_length as the number of compressed bytes
        compressed_bytes = read(io, asset.encoded_length)
        
        println("DEBUG: Decompressing data - offset=$(asset.offset), compressed_bytes=$(length(compressed_bytes))")
        
        local decompressed_bytes
        try
            decompressed_bytes = Libz.inflate(compressed_bytes)
            println("DEBUG: Decompression successful - decompressed_bytes=$(length(decompressed_bytes))")
        catch e
            @error "ZLIB DECOMPRESSION FAILED. This is likely due to an incorrect offset or corrupt data in the .ibd file."
            @error "Asset offset: $(asset.offset), Encoded length: $(asset.encoded_length)"
            # Print first 16 bytes to stderr for diagnosis
            bytes_to_print = min(16, length(compressed_bytes))
            @error "First $bytes_to_print bytes of the data chunk we tried to decompress:"
            println(stderr, view(compressed_bytes, 1:bytes_to_print))
            rethrow(e)
        end
        
        # Use an IOBuffer to safely read the data
        bytes_io = IOBuffer(decompressed_bytes)
        n_elements = bytes_io.size ÷ sizeof(format)
        array = Array{format}(undef, n_elements)
        read!(bytes_io, array)
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
function _iterate_uncompressed_fast(f::Function, data::MSIData, source::ImzMLSource)
    # Optimized path for uncompressed data using buffer reuse
    max_points = maximum(meta -> meta.mz_asset.encoded_length, data.spectra_metadata)
    mz_buffer = Vector{source.mz_format}(undef, max_points)
    int_buffer = Vector{source.intensity_format}(undef, max_points)

    for i in 1:length(data.spectra_metadata)
        meta = data.spectra_metadata[i]
        nPoints = meta.mz_asset.encoded_length

        if nPoints == 0
            f(i, view(mz_buffer, 0:-1), view(int_buffer, 0:-1))
            continue
        end

        mz_view = view(mz_buffer, 1:nPoints)
        int_view = view(int_buffer, 1:nPoints)

        if meta.mz_asset.offset < meta.int_asset.offset
            seek(source.ibd_handle, meta.mz_asset.offset)
            read!(source.ibd_handle, mz_view)
            read!(source.ibd_handle, int_view)
        else
            seek(source.ibd_handle, meta.int_asset.offset)
            read!(source.ibd_handle, int_view)
            read!(source.ibd_handle, mz_view)
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
function _iterate_compressed_fast(f::Function, data::MSIData, source::ImzMLSource)
    # Path for datasets containing at least one compressed spectrum.
    # This path reads and decompresses each spectrum individually.
    for i in 1:length(data.spectra_metadata)
        meta = data.spectra_metadata[i]
        
        if meta.mz_asset.encoded_length == 0 && meta.int_asset.encoded_length == 0
            f(i, source.mz_format[], source.intensity_format[])
            continue
        end

        # Read and decompress each array
        mz_array = read_compressed_array(source.ibd_handle, meta.mz_asset, source.mz_format)
        intensity_array = read_compressed_array(source.ibd_handle, meta.int_asset, source.intensity_format)
        
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
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource)
    if isempty(data.spectra_metadata)
        return
    end
    
    # Check if ANY spectra are compressed and dispatch to the appropriate implementation
    any_compressed = any(meta -> meta.mz_asset.is_compressed || meta.int_asset.is_compressed, 
                        data.spectra_metadata)
    
    if any_compressed
        _iterate_compressed_fast(f, data, source)
    else
        _iterate_uncompressed_fast(f, data, source)
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