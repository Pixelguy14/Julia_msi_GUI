# src/MSIData.jl

"""
This file defines the unified MSIData object and the associated data access layer,
including caching and iteration logic, for handling large mzML and imzML datasets
efficiently.
"""

using Base64, Libz, Serialization # For reading binary data

# Abstract type for different data sources (e.g., mzML, imzML)
# This allows dispatching to the correct binary reading logic.
abstract type MSDataSource end

# Concrete type for imzML data sources
struct ImzMLSource <: MSDataSource
    ibd_handle::IO
    mz_format::Type
    intensity_format::Type
end

# Concrete type for mzML data sources
struct MzMLSource <: MSDataSource
    file_handle::IO
    mz_format::Type
    intensity_format::Type
end

# Struct to hold info about a binary data array (mz or intensity)
struct SpectrumAsset
    format::Type
    is_compressed::Bool
    offset::Int64
    encoded_length::Int32
    # For mzML, axis_type is needed to distinguish mz from intensity.
    # For imzML, this can be ignored as the order is fixed.
    axis_type::Symbol
end

# Metadata for a single spectrum, common to both formats
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

# The main unified data object
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
function read_spectrum_from_disk(source::ImzMLSource, meta::SpectrumMetadata)
    # For imzML, the binary data is raw, not base64 encoded. 
    # The `encoded_length` field in this case holds the number of points.
    mz = Array{source.mz_format}(undef, meta.mz_asset.encoded_length)
    intensity = Array{source.intensity_format}(undef, meta.int_asset.encoded_length)

    seek(source.ibd_handle, meta.mz_asset.offset)
    read!(source.ibd_handle, mz)

    seek(source.ibd_handle, meta.int_asset.offset)
    read!(source.ibd_handle, intensity)

    return mz, intensity
end

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
    println("Calculating total spectrum (2-pass method for imzML)...")

    # 1. First Pass: Find the global m/z range
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

    if !isfinite(global_min_mz)
        @warn "Could not determine a valid m/z range. All spectra might be empty."
        return (Float64[], Float64[])
    end
    println("  Global m/z range found: [$(global_min_mz), $(global_max_mz)]")

    # 2. Define Bins
    mz_bins = range(global_min_mz, stop=global_max_mz, length=num_bins)
    intensity_sum = zeros(Float64, num_bins)
    bin_step = step(mz_bins)

    # 3. Second Pass: Sum intensities into bins
    println("  Pass 2: Summing intensities into $num_bins bins...")
    _iterate_spectra_fast(msi_data) do idx, mz, intensity
        if isempty(mz)
            return # continue equivalent
        end
        for i in eachindex(mz)
            bin_index = clamp(round(Int, (mz[i] - global_min_mz) / bin_step) + 1, 1, num_bins)
            intensity_sum[bin_index] += intensity[i]
        end
    end
    
    println("Total spectrum calculation complete.")
    return (collect(mz_bins), intensity_sum)
end

function get_total_spectrum_mzml(msi_data::MSIData; num_bins::Int=2000)
    println("Calculating total spectrum (single-pass optimization for mzML)...")
    num_spectra = length(msi_data.spectra_metadata)
    if num_spectra == 0
        return (Float64[], Float64[])
    end

    temp_path, temp_io = mktemp()
    try
        # --- Pass 1: Write spectra to temp file and find min/max m/z ---
        println("  Pass 1: Caching spectra and finding global m/z range...")
        global_min_mz = Inf
        global_max_mz = -Inf
        
        _iterate_spectra_fast(msi_data) do idx, mz, intensity
            Serialization.serialize(temp_io, (mz, intensity))
            if !isempty(mz)
                local_min, local_max = extrema(mz)
                global_min_mz = min(global_min_mz, local_min)
                global_max_mz = max(global_max_mz, local_max)
            end
        end

        flush(temp_io)

        if !isfinite(global_min_mz)
            @warn "Could not determine a valid m/z range. All spectra might be empty."
            return (Float64[], Float64[])
        end
        println("  Global m/z range found: [$(global_min_mz), $(global_max_mz)]")

        # --- Pass 2: Read from temp file and bin intensities ---
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
        
        println("Total spectrum calculation complete.")
        return (collect(mz_bins), intensity_sum)

    finally
        close(temp_io)
        rm(temp_path, force=true)
    end
end

"""
    get_total_spectrum(msi_data::MSIData; num_bins::Int=20000) -> Tuple{Vector{Float64}, Vector{Float64}}

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

struct MSIDataIterator
    data::MSIData
end

"""
    IterateSpectra(data::MSIData)

Returns an iterator that yields each spectrum, processing the file sequentially
with minimal memory overhead.

This function is the core of the "Event-driven" access pattern.
"""
function IterateSpectra(data::MSIData)
    return MSIDataIterator(data)
end

# Define the iteration interface for the custom iterator
Base.length(it::MSIDataIterator) = length(it.data.spectra_metadata)
Base.eltype(it::MSIDataIterator) = Tuple{Int, Tuple{Vector, Vector}} # Index, (mz, intensity)

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
    _iterate_spectra_fast(f::Function, data::MSIData)

An internal, high-performance iterator for bulk processing.
It avoids the overhead of `GetSpectrum` by reading directly from the file stream
and reusing pre-allocated buffers.

- `f`: A function to call for each spectrum, with signature `f(index, mz, intensity)`.
- `data`: The `MSIData` object.
"""
function _iterate_spectra_fast_impl(f::Function, data::MSIData, source::ImzMLSource)
    # This implementation is for imzML and is optimized for performance by
    # using pre-allocated buffers.
    max_points = 0
    for meta in data.spectra_metadata
        # For imzML, encoded_length is the number of points
        max_points = max(max_points, meta.mz_asset.encoded_length)
    end

    # If all spectra are empty, just call f with empty arrays
    if max_points == 0 && !isempty(data.spectra_metadata)
        for i in 1:length(data.spectra_metadata)
            f(i, source.mz_format[], source.intensity_format[])
        end
        return
    end

    mz_buffer = Vector{source.mz_format}(undef, max_points)
    int_buffer = Vector{source.intensity_format}(undef, max_points)

    for i in 1:length(data.spectra_metadata)
        meta = data.spectra_metadata[i]
        nPoints = meta.mz_asset.encoded_length

        if nPoints == 0
            f(i, source.mz_format[], source.intensity_format[])
            continue
        end

        mz_view = view(mz_buffer, 1:nPoints)
        int_view = view(int_buffer, 1:nPoints)

        seek(source.ibd_handle, meta.mz_asset.offset)
        read!(source.ibd_handle, mz_view)

        seek(source.ibd_handle, meta.int_asset.offset)
        read!(source.ibd_handle, int_view)

        f(i, mz_view, int_view)
    end
end

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

function _iterate_spectra_fast(f::Function, data::MSIData)
    # Dispatch to the correct implementation based on the source type
    _iterate_spectra_fast_impl(f, data, data.source)
end
