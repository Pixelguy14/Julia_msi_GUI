# /home/pixel/Documents/Cinvestav_2025/JuliaMSI/src/mzML.jl

# This file is responsible for parsing metadata from .mzML files.
# It has been refactored to produce a unified MSIData object.

# This function is kept internal to the mzML parsing process
function get_spectrum_asset_metadata(stream)
    start_pos = position(stream)
    
    bda_tag = find_tag(stream, r"<binaryDataArray\s+encodedLength=\"(\d+)\"" )
    if bda_tag === nothing
        error("Cannot find binaryDataArray")
    end
    encoded_length = parse(Int32, bda_tag.captures[1])

    # Initialize parameters as separate variables
    data_format = Float64
    compression_flag = false
    axis = :mz

    while true
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
            # println("DEBUG: Processing accession: $acc_str")
            
            # Direct inline logic - NO FUNCTION CALLS
            if acc_str == "MS:1000514"
                axis = :mz
                # println("DEBUG: Set axis_type to :mz")
            elseif acc_str == "MS:1000515"
                axis = :intensity
                # println("DEBUG: Set axis_type to :intensity")
            elseif acc_str == "MS:1000518"
                data_format = Int16
                # println("DEBUG: Set format to Int16")
            elseif acc_str == "MS:1000519"
                data_format = Int32
                # println("DEBUG: Set format to Int32")
            elseif acc_str == "MS:1000520"
                data_format = Float64
                # println("DEBUG: Set format to Float64")
            elseif acc_str == "MS:1000521"
                data_format = Float32
                # println("DEBUG: Set format to Float32")
            elseif acc_str == "MS:1000522"
                data_format = Int64
                # println("DEBUG: Set format to Int64")
            elseif acc_str == "MS:1000523"
                data_format = Float64
                # println("DEBUG: Set format to Float64")
            elseif acc_str == "MS:1000574"
                compression_flag = true
                # println("DEBUG: Set is_compressed to true")
            elseif acc_str == "MS:1000576"
                compression_flag = false
                # println("DEBUG: Set is_compressed to false")
            else
                # println("DEBUG: Unknown accession: $acc_str")
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
function parse_spectrum_metadata(stream, offset::Int64)
    seek(stream, offset)
    
    id_match = find_tag(stream, r"<spectrum\s+index=\"\d+\"\s+id=\"([^\"]+)" )
    id = id_match === nothing ? "" : id_match.captures[1]

    asset1 = get_spectrum_asset_metadata(stream)
    asset2 = get_spectrum_asset_metadata(stream)

    mz_asset = asset1.axis_type == :mz ? asset1 : asset2
    int_asset = asset1.axis_type == :intensity ? asset1 : asset2

    # Create the new unified metadata object
    # For mzML, x and y coordinates are not applicable, so we use 0.
    return SpectrumMetadata(0, 0, id, mz_asset, int_asset)
end

function parse_offset_list(stream)
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

# This is the main lazy-loading function for mzML, now returning an MSIData object.
function load_mzml_lazy(file_path::String; cache_size::Int=100)
    println("DEBUG: Opening file stream for $file_path")
    stream = open(file_path, "r")

    try
        println("DEBUG: Seeking to end of file.")
        seekend(stream)
        println("DEBUG: Skipping back to read footer.")
        skip(stream, -4096)
        footer = read(stream, String)
        println("DEBUG: Footer read successfully.")
        
        index_offset_match = match(r"<indexListOffset>(\d+)</indexListOffset>", footer)
        if index_offset_match === nothing
            close(stream)
            error("Could not find <indexListOffset>. File may not be an indexed mzML.")
        end
        println("DEBUG: Found indexListOffset.")
        
        index_offset = parse(Int64, index_offset_match.captures[1])
        println("DEBUG: Seeking to index list at offset $index_offset.")
        seek(stream, index_offset)
        
        println("DEBUG: Searching for '<index name=\"spectrum\">'.")
        if find_tag(stream, r"<index\s+name=\"spectrum\"") === nothing
            close(stream)
            error("Could not find spectrum index.")
        end
        println("DEBUG: Found spectrum index tag. ")

        println("DEBUG: Parsing spectrum offsets...")
        spectrum_offsets = parse_offset_list(stream)
        if isempty(spectrum_offsets)
            close(stream)
            error("No spectrum offsets found.")
        end
        println("DEBUG: Found $(length(spectrum_offsets)) spectrum offsets.")

        println("DEBUG: Parsing metadata for each spectrum...")
        spectra_metadata = [parse_spectrum_metadata(stream, offset) for offset in spectrum_offsets]
        println("DEBUG: Metadata parsing complete.")

        # Assuming uniform data formats, take from the first spectrum
        first_meta = spectra_metadata[1]
        mz_format = first_meta.mz_asset.format
        intensity_format = first_meta.int_asset.format
        
        source = MzMLSource(stream, mz_format, intensity_format)
        println("DEBUG: Creating MSIData object.")
        return MSIData(source, spectra_metadata, (0, 0), nothing, cache_size)

    catch e
        close(stream) # Ensure stream is closed on error
        rethrow(e)
    end
end

"""
    LoadMzml(fileName::String)

Eagerly loads all spectra from a .mzML file.
Provided for backward compatibility. Now uses the new MSIData architecture.
"""
function LoadMzml(fileName::String)
    # Use the lazy loader to get the MSIData object
    msi_data = load_mzml_lazy(fileName, cache_size=0) # No need to cache if we load all
    
    num_spectra = length(msi_data.spectra_metadata)
    spectra_matrix = Array{Any}(undef, (2, num_spectra))
    
    for i in 1:num_spectra
        # Use the new GetSpectrum API
        mz, intensity = GetSpectrum(msi_data, i)
        spectra_matrix[1, i] = mz
        spectra_matrix[2, i] = intensity
    end
    
    # The finalizer on MSIData will close the handle, so no need to close it here.
    
    return spectra_matrix
end
