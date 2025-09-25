# /home/pixel/Documents/Cinvestav_2025/JuliaMSI/src/mzML_lazy.jl

using Base64, Libz
include("ParserHelpers.jl")

# ============================================================================
# Data Structures for Lazy Loading
# ============================================================================ 

struct SpectrumAsset
    format::Type
    is_compressed::Bool
    offset::Int64
    encoded_length::Int32
    axis_type::Symbol # :mz or :intensity
end

struct MzMLSpectrumMetadata
    id::String
    mz_asset::SpectrumAsset
    intensity_asset::SpectrumAsset
end

mutable struct MzMLData
    file_handle::IO
    spectra_metadata::Vector{MzMLSpectrumMetadata}

    function MzMLData(file_handle::IO, spectra_metadata::Vector{MzMLSpectrumMetadata})
        obj = new(file_handle, spectra_metadata)
        finalizer(obj) do o
            if isopen(o.file_handle)
                close(o.file_handle)
            end
        end
        return obj
    end
end

# ============================================================================ 
# mzML Parser Implementation
# ============================================================================ 

function get_spectrum_asset_metadata(stream)
    start_pos = position(stream)
    
    bda_tag = find_tag(stream, r"<binaryDataArray\s+encodedLength=\"(\d+)\"")
    if bda_tag === nothing; error("Cannot find binaryDataArray"); end
    encoded_length = parse(Int32, bda_tag.captures[1])

    params = CVParams(Float64, false, :mz)

    while true
        line = readline(stream)
        if occursin("</binaryDataArray>", line); break; end
        if occursin("<cvParam", line)
            accession = get_attribute(line, "accession")
            if accession === nothing; continue; end
            update_cv_params!(params, accession.captures[1])
        end
    end

    seek(stream, start_pos)
    binary_tag = find_tag(stream, r"<binary>")
    if binary_tag === nothing; error("Cannot find binary tag"); end
    binary_offset = position(stream)

    find_tag(stream, r"</binaryDataArray>") # Position for next call

    return SpectrumAsset(params.format, params.is_compressed, binary_offset, encoded_length, params.axis_type)
end

function parse_spectrum_metadata(stream, offset::Int64)
    seek(stream, offset)
    
    id_match = find_tag(stream, r"<spectrum\s+index=\"\d+\"\s+id=\"([^"]+)\"")
    id = id_match === nothing ? "" : id_match.captures[1]

    asset1 = get_spectrum_asset_metadata(stream)
    asset2 = get_spectrum_asset_metadata(stream)

    mz_asset = asset1.axis_type == :mz ? asset1 : asset2
    int_asset = asset1.axis_type == :intensity ? asset1 : asset2

    return MzMLSpectrumMetadata(id, mz_asset, int_asset)
end

function parse_offset_list(stream)
    offsets = Int64[]
    while true
        line = readline(stream)
        m = match(r"<offset\s+idRef=\"\S+\">(\d+)</offset>", line)
        if m !== nothing
            push!(offsets, parse(Int64, m.captures[1]))
        elseif occursin("</index>", line) || occursin("</indexedmzML>", line)
            break
        end
    end
    return offsets
end

function load_mzml_lazy(file_path::String)
    stream = open(file_path)

    try
        seekend(stream)
        skip(stream, -4096) # Read last 4KB
        footer = read(stream, String)
        
        index_offset_match = match(r"<indexListOffset>(\d+)</indexListOffset>", footer)
        if index_offset_match === nothing
            error("Could not find <indexListOffset>. File may not be an indexed mzML.")
        end
        
        index_offset = parse(Int64, index_offset_match.captures[1])
        seek(stream, index_offset)
        
        if find_tag(stream, r"<index\s+name=\"spectrum\"") === nothing
            error("Could not find spectrum index.")
        end

        spectrum_offsets = parse_offset_list(stream)
        if isempty(spectrum_offsets)
            error("No spectrum offsets found.")
        end

        spectra_metadata = [parse_spectrum_metadata(stream, offset) for offset in spectrum_offsets]
        
        # We return the stream inside MzMLData, so we don't close it here.
        # The finalizer on MzMLData will handle it.
        return MzMLData(stream, spectra_metadata)

    catch e
        close(stream) # Close file on error
        rethrow(e)
    end
end

# ============================================================================ 
# Data Access Functions
# ============================================================================ 

function read_binary_vector(io::IO, asset::SpectrumAsset)
    seek(io, asset.offset)
    
    # Read the raw base64 data. The regex is slow. Reading bytes is better.
    # The old code did `Vector{UInt8}(base64Vec.captures[1])`. This implies ASCII.
    raw_b64 = read(io, asset.encoded_length)
    
    decoded_bytes = base64decode(String(raw_b64))
    
    bytes_io = IOBuffer(asset.is_compressed ? Libz.inflate(decoded_bytes) : decoded_bytes)
    
    n_elements = Int(bytes_io.size / sizeof(asset.format))
    out_array = Array{asset.format}(undef, n_elements)
    read!(bytes_io, out_array)
    
    return out_array
end

function get_spectrum(data::MzMLData, index::Int)
    meta = data.spectra_metadata[index]
    
    mz_array = read_binary_vector(data.file_handle, meta.mz_asset)
    intensity_array = read_binary_vector(data.file_handle, meta.intensity_asset)
    
    return (id=meta.id, mz=mz_array, intensity=intensity_array)
end

"""
    LoadMzml(fileName::String)

Eagerly loads all spectra from a .mzML file.
Provided for backward compatibility.
"""
function LoadMzml(fileName::String)
    mzml_data = load_mzml_lazy(fileName)
    
    num_spectra = length(mzml_data.spectra_metadata)
    spectra_matrix = Array{Any}(undef, (2, num_spectra))
    
    for i in 1:num_spectra
        spectrum = get_spectrum(mzml_data, i)
        spectra_matrix[1, i] = spectrum.mz
        spectra_matrix[2, i] = spectrum.intensity
    end
    
    # Eager load, so we can close the handle now.
    close(mzml_data.file_handle)
    
    return spectra_matrix
end
