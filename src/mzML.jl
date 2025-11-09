# src/mzML.jl

# This file is responsible for parsing metadata from .mzML files.
# It has been refactored to produce a unified MSIData object.

# Constants for CV parameter accessions - defined once for performance
const MZ_AXIS_ACCESSION = "MS:1000514"
const INTENSITY_AXIS_ACCESSION = "MS:1000515"
const COMPRESSION_ACCESSION = "MS:1000574"
const NO_COMPRESSION_ACCESSION = "MS:1000576"

# Data format accessions as constants
const DATA_FORMAT_ACCESSIONS = Dict{String, DataType}(
    "MS:1000518" => Int16,
    "MS:1000519" => Int32,
    "MS:1000520" => Float64,
    "MS:1000521" => Float32,
    "MS:1000522" => Int64,
    "MS:1000523" => Float64
)

"""
    get_spectrum_asset_metadata(stream::IO)

Parses a `<binaryDataArray>` block within an mzML file to extract metadata
for a single data array (e.g., m/z or intensity). It reads CV parameters to
determine the data type, compression, and axis type.

# Arguments
- `stream`: An IO stream positioned at the beginning of a `<binaryDataArray>` block.

# Returns
- A `SpectrumAsset` struct containing the parsed metadata, including the binary
  data offset, encoded length, format, and compression status.
"""
function get_spectrum_asset_metadata(stream::IO)
    start_pos = position(stream)
    
        bda_tag = find_tag(stream, r"<binaryDataArray\s+encodedLength=\"(\d+)\"")
    
        if bda_tag === nothing
    
            throw(FileFormatError("Cannot find binaryDataArray"))
    
        end
    encoded_length = parse(Int32, bda_tag.captures[1])

    # Initialize parameters as separate variables with concrete types
    data_format::DataType = Float64
    compression_flag::Bool = false
    axis::Symbol = :mz

    while !eof(stream)
        line = readline(stream)
        if occursin("</binaryDataArray>", line)
            break
        end
        if occursin("<cvParam", line)
            accession = get_attribute(line, "accession")
            if accession === nothing
                continue
            end
            acc_str = accession.captures[1]
            
            # Use constant comparisons and dictionary lookup for better performance
            if acc_str == MZ_AXIS_ACCESSION
                axis = :mz
            elseif acc_str == INTENSITY_AXIS_ACCESSION
                axis = :intensity
            elseif haskey(DATA_FORMAT_ACCESSIONS, acc_str)
                data_format = DATA_FORMAT_ACCESSIONS[acc_str]
            elseif acc_str == COMPRESSION_ACCESSION
                compression_flag = true
            elseif acc_str == NO_COMPRESSION_ACCESSION
                compression_flag = false
            end
        end
    end

    seek(stream, start_pos)
    readuntil(stream, "<binary>")
    binary_offset = position(stream)

    # Move stream to the end of the binary data array for the next iteration
    readuntil(stream, "</binaryDataArray>")

    # Create SpectrumAsset directly from the variables
    return SpectrumAsset(data_format, compression_flag, binary_offset, encoded_length, axis)
end

# This function is updated to return the generic SpectrumMetadata struct
"""
    parse_spectrum_metadata(stream::IO, offset::Int64)

Parses an entire `<spectrum>` block from an mzML file, given a starting offset.
It extracts the spectrum ID and calls `get_spectrum_asset_metadata` to parse
the m/z and intensity array metadata.

# Arguments
- `stream`: An IO stream for the mzML file.
- `offset`: The byte offset where the `<spectrum>` block begins.

# Returns
- A `SpectrumMetadata` struct containing the parsed metadata for one spectrum.
"""
function parse_spectrum_metadata(stream::IO, offset::Int64)
    seek(stream, offset)
    
    # Read the whole spectrum block to parse mode
    spectrum_start_pos = position(stream)
    line = ""
    spectrum_buffer = IOBuffer()
    while !eof(stream)
        line = readline(stream)
        write(spectrum_buffer, line)
        if occursin("</spectrum>", line)
            break
        end
    end
    spectrum_xml = String(take!(spectrum_buffer))
    seek(stream, spectrum_start_pos) # Reset for other parsing

    id_match = match(r"<spectrum\s+index=\"\d+\"\s+id=\"([^\"]+)", spectrum_xml)
    id = id_match === nothing ? "" : id_match.captures[1]

    # Determine mode from the XML block
    mode = UNKNOWN
    if occursin("MS:1000127", spectrum_xml)
        mode = CENTROID
    elseif occursin("MS:1000128", spectrum_xml)
        mode = PROFILE
    end

    # Find where the binary data list starts to parse assets
    binary_list_match = findfirst("<binaryDataArrayList", spectrum_xml)
    if binary_list_match !== nothing
        seek(stream, spectrum_start_pos + binary_list_match.start - 1)
    end

    asset1 = get_spectrum_asset_metadata(stream)
    asset2 = get_spectrum_asset_metadata(stream)

    # Determine which asset is m/z and which is intensity
    mz_asset, int_asset = if asset1.axis_type == :mz
        (asset1, asset2)
    else
        (asset2, asset1)
    end

    # Create the new unified metadata object
    # For mzML, x and y coordinates are not applicable, so we use 0.
    return SpectrumMetadata(Int32(0), Int32(0), id, mode, mz_asset, int_asset)
end

"""
    parse_offset_list(stream::IO)

Parses the `<index name="spectrum">` block in an indexed mzML file to extract
the byte offsets for each spectrum.

# Arguments
- `stream`: An IO stream positioned at the start of the `<index>` block.

# Returns
- A `Vector{Int64}` containing the byte offsets for all spectra.
"""
function parse_offset_list(stream::IO)
    offsets = Int64[]
    offset_regex = r"<offset[^>]*>(\d+)</offset>"
    
    # Actively search for offset tags, ignoring other lines until the end of the index is found.
    while !eof(stream)
        line = readline(stream)
        
        # First, check for the end condition.
        if occursin("</index>", line) || occursin("</indexedmzML>", line)
            break
        end

        # If it's not the end, see if it's an offset tag.
        m = match(offset_regex, line)
        if m !== nothing
            push!(offsets, parse(Int64, m.captures[1]))
        end
        # If it's neither, ignore the line and continue the loop.
    end
    
    return offsets
end

"""
    find_index_offset(stream::IO)::Int64

Finds the index offset in an mzML file by reading from the end.
Optimized version with better buffer management.
"""
function find_index_offset(stream::IO)::Int64
    file_size = filesize(stream)
    seekend(stream)
    
    # Read larger chunk for better chance of finding the offset
    chunk_size = min(8192, file_size)
    seek(stream, file_size - chunk_size)
    footer = read(stream, String)
    
    index_offset_match = match(r"<indexListOffset>(\d+)</indexListOffset>", footer)
    if index_offset_match === nothing
        throw(FileFormatError("Could not find <indexListOffset>. File may not be an indexed mzML."))
    end
    
    return parse(Int64, index_offset_match.captures[1])
end

# This is the main lazy-loading function for mzML, now returning an MSIData object.
"""
    load_mzml_lazy(file_path::String; cache_size::Int=100)

Lazily loads an indexed `.mzML` file by parsing only the metadata. It reads the
spectrum index from the end of the file to get the offsets of each spectrum,
then parses the metadata for each spectrum without loading the binary data.

# Arguments
- `file_path`: The path to the `.mzML` file.
- `cache_size`: The number of spectra to hold in an LRU cache for faster access.

# Returns
- An `MSIData` object ready for lazy data access.
"""
function load_mzml_lazy(file_path::String; cache_size::Int=100)
    println("DEBUG: Opening file stream for $file_path")
    ts_stream = ThreadSafeFileHandle(file_path, "r")

    try
        println("DEBUG: Finding index offset...")
        index_offset = find_index_offset(ts_stream.handle)
        println("DEBUG: Seeking to index list at offset $index_offset.")
        seek(ts_stream.handle, index_offset)
        
        println("DEBUG: Searching for '<index name=\"spectrum\">'.")
        if find_tag(ts_stream.handle, r"<index\s+name=\"spectrum\"") === nothing
            throw(FileFormatError("Could not find spectrum index."))
        end
        println("DEBUG: Found spectrum index tag.")

        println("DEBUG: Parsing spectrum offsets...")
        spectrum_offsets = parse_offset_list(ts_stream.handle)
        if isempty(spectrum_offsets)
            throw(FileFormatError("No spectrum offsets found."))
        end
        num_spectra = length(spectrum_offsets)
        println("DEBUG: Found $num_spectra spectrum offsets.")

        println("DEBUG: Parsing metadata for each spectrum...")
        # Pre-allocate the metadata vector for better performance
        spectra_metadata = Vector{SpectrumMetadata}(undef, num_spectra)
        
        # Use @inbounds for faster indexing in the loop
        @inbounds for i in 1:num_spectra
            spectra_metadata[i] = parse_spectrum_metadata(ts_stream.handle, spectrum_offsets[i])
            
            # Progress reporting for large files
            if i % 1000 == 0
                println("DEBUG: Processed $i/$num_spectra spectra")
            end
        end
        println("DEBUG: Metadata parsing complete.")

        # Assuming uniform data formats, take from the first spectrum
        first_meta = spectra_metadata[1]
        mz_format = first_meta.mz_asset.format
        intensity_format = first_meta.int_asset.format
        
        source = MzMLSource(ts_stream, mz_format, intensity_format)
        println("DEBUG: Creating MSIData object.")
        return MSIData(source, spectra_metadata, (0, 0), nothing, cache_size)

    catch e
        close(ts_stream) # Ensure stream is closed on error
        rethrow(e)
    end
end

#=
"""
    LoadMzml(fileName::String)

Eagerly loads all spectra from a .mzML file into memory.
This function now uses the new lazy-loading
architecture internally but presents the data in the old format.

# Arguments
- `fileName`: The path to the `.mzML` file.

# Returns
- A `2xN` matrix where `N` is the number of spectra. The first row contains
  m/z arrays and the second row contains intensity arrays.
"""
function LoadMzml(fileName::String)
    # Use the lazy loader to get the MSIData object
    msi_data = load_mzml_lazy(fileName, cache_size=0) # No need to cache if we load all
    
    num_spectra = length(msi_data.spectra_metadata)
    
    # FIXED: Use concrete typed arrays instead of Array{Any}
    spectra_matrix = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, num_spectra)
    
    # Pre-allocate and use bounds checking optimization
    @inbounds for i in 1:num_spectra
        # Use the new GetSpectrum API
        mz, intensity = GetSpectrum(msi_data, i)
        spectra_matrix[i] = (mz, intensity)
    end
    
    # Convert to the expected 2xN format if needed by downstream code
    # Note: This maintains the original interface but with better typing
    result_matrix = Array{Any}(undef, (2, num_spectra))
    @inbounds for i in 1:num_spectra
        result_matrix[1, i] = spectra_matrix[i][1]
        result_matrix[2, i] = spectra_matrix[i][2]
    end
    
    return result_matrix
end
=#

#=
"""
    load_mzml_batch(file_path::String, spectrum_indices::AbstractVector{Int})

Loads a specific batch of spectra from an mzML file efficiently.
Useful for parallel processing or when only specific spectra are needed.

# Arguments
- `file_path`: Path to the mzML file
- `spectrum_indices`: Indices of spectra to load

# Returns
- Vector of (mz_array, intensity_array) tuples
"""
function load_mzml_batch(file_path::String, spectrum_indices::AbstractVector{Int})
    msi_data = load_mzml_lazy(file_path)
    num_to_load = length(spectrum_indices)
    
    # Pre-allocate result with concrete types
    results = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, num_to_load)
    
    @inbounds for (i, idx) in enumerate(spectrum_indices)
        mz, intensity = GetSpectrum(msi_data, idx)
        results[i] = (mz, intensity)
    end
    
    return results
end
=#

#=
"""
    get_mzml_summary(file_path::String)::NamedTuple

Quickly extracts summary information from an mzML file without loading all metadata.

# Returns
- Named tuple with: num_spectra, mz_range, intensity_range, data_formats
"""
function get_mzml_summary(file_path::String)::NamedTuple
    ts_stream = ThreadSafeFileHandle(file_path, "r")
    try
        index_offset = find_index_offset(ts_stream.handle)
        seek(ts_stream.handle, index_offset)
        
        if find_tag(ts_stream.handle, r"<index\s+name=\"spectrum\"") === nothing
            throw(FileFormatError("Could not find spectrum index."))
        end
        
        spectrum_offsets = parse_offset_list(ts_stream.handle)
        num_spectra = length(spectrum_offsets)
        
        # Sample first few spectra to determine formats and ranges
        sample_size = min(10, num_spectra)
        sample_metadata = [parse_spectrum_metadata(ts_stream.handle, spectrum_offsets[i]) for i in 1:sample_size]
        
        # Extract formats from sample
        mz_format = sample_metadata[1].mz_asset.format
        intensity_format = sample_metadata[1].int_asset.format
        
        return (num_spectra=num_spectra, 
                mz_format=mz_format, 
                intensity_format=intensity_format,
                sample_spectra=sample_size)
    finally
        close(ts_stream)
    end
end
=#

#=
# Performance optimization: Cache frequently used regex patterns
const PRE_COMPILED_REGEX = (
    encoded_length = r"<binaryDataArray\s+encodedLength=\"(\d+)\"",
    spectrum_id = r"<spectrum\s+index=\"\d+\"\s+id=\"([^\"]+)",
    offset = r"<offset[^>]*>(\d+)</offset>",
    index_list = r"<index\s+name=\"spectrum\"",
    index_offset = r"<indexListOffset>(\d+)</indexListOffset>"
)
=#
