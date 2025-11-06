# src/imzML.jl
using Images, Statistics, CairoMakie, DataFrames, Printf, ColorSchemes, StatsBase

"""
This file provides a library for parsing `.imzML` and `.ibd` files in pure Julia.

Core Functions:
- `load_imzml_lazy`: The main function that orchestrates the parsing.
- Helper functions for reading XML metadata and binary spectral data.
"""


# ============================================================================
#
#
# imzML Parser Implementation
#
# ============================================================================

"""
    axes_config_img(stream)

Determines the storage order of the m/z and intensity arrays.
"""
function axes_config_img(stream::IO)
    param_groups = Dict{String, SpecDim}()
    find_tag(stream, r"<referenceableParamGroupList")

    while true
        pos = position(stream)
        line = readline(stream)

        if eof(stream) || occursin("</referenceableParamGroupList>", line)
            break
        end

        id_match = match(r"<referenceableParamGroup id=\"([^\"]+)\"", line)
        if id_match !== nothing
            id = id_match.captures[1]
            spec_dim = configure_spec_dim(stream)
            param_groups[id] = spec_dim
        end
    end
    return param_groups
end

"""
    get_img_dimensions(stream)

Reads the maximum X and Y dimensions and total spectrum count from the metadata.
"""
function get_img_dimensions(stream::IO)
    find_tag(stream, r"^\s*<(scanSettings )")
    n = 2
    dim = zeros(Int32, 3)

    while n > 0 && !eof(stream)
        currLine = readline(stream)
        if occursin("<cvParam", currLine)
            accession_match = get_attribute(currLine, "accession")
            if accession_match !== nothing
                accession = accession_match.captures[1]
                if accession == "IMS:1000042" || accession == "IMS:1000043"
                    value_match = get_attribute(currLine, "value")
                    if value_match !== nothing
                        axis_idx = (accession == "IMS:1000042") ? 1 : 2
                        dim[axis_idx] = parse(Int32, value_match.captures[1])
                        n -= 1
                    end
                end
            end
        end
    end

    count_tag = find_tag(stream, r"^\s*<spectrumList(.+)")
    count_match = get_attribute(count_tag.captures[1], "count")
    dim[3] = parse(Int32, count_match.captures[1])
    return dim
end

"""
    get_spectrum_tag_offset(stream)

Calculates the character offset within a `<spectrum>` tag, ignoring attribute values.
"""
function get_spectrum_tag_offset(stream::IO)
    offset = position(stream)
    tag = find_tag(stream, r"^\s*<spectrum (.+)")
    first = 1

    while true
        value = match(r"[^=]+\"([^\"]+)\"", tag.captures[1][first:end])
        if value === nothing
            break
        end
        first += value.offsets[1] + length(value.captures[1])
        offset += length(value.captures[1])
    end
    return offset
end

"""
    get_spectrum_attributes(stream, hIbd)

Reads metadata to determine the byte offsets and data types for reading spectra.
"""
function get_spectrum_attributes(stream::IO, hIbd::IO)
    skip = Vector{UInt32}(undef, 8)
    offset = get_spectrum_tag_offset(stream)

    tag = find_tag(stream, r" accession=\"IMS:100005(\d)\"(.+)")
    skip[1] = tag.captures[1][1] - '0' + 1
    skip[2] = xor(skip[1], 3)

    value = match(r" value=\"(\d+)\".+", tag.captures[2])
    skip[5] = position(stream) - offset - length(value.match) - 2
    value = match(r" value=\"\d+\".+", readline(stream))
    skip[6] = value.offset

    offset = position(stream)
    tag = find_tag(stream, r"^\s*<referenceableParamGroupRef(.+)")
    value = get_attribute(tag.captures[1], "ref")
    skip[3] = (value.captures[1] == "intensityArray") + 3
    skip[4] = xor(skip[3], 7)

    k = 2
    while k != 0
        tag = find_tag(stream, r" accession=\"IMS:100010(\d)\"(.+)")
        accession_val = tag.captures[1][1]
        if accession_val == '2'
            value = get_attribute(tag.match, "value")
            seek(hIbd, parse(Int64, value.captures[1]))
            k -= 1
        elseif accession_val == '3'
            skip[7] = position(stream) - offset - length(tag.match)
            offset = position(stream)
            k -= 1
        end
    end

    find_tag(stream, r"^\s*</spectrum>")
    skip[8] = position(stream) - offset
    return skip
end

function get_spectrum_attributes_v2(stream::IO, hIbd::IO)
    # Look for position x
    readuntil(stream, "IMS:1000050")
    x_skip = 0
    
    # Look for position y  
    readuntil(stream, "IMS:1000051")
    y_skip = 0
    
    # Determine order of mz/intensity arrays
    pos_before = position(stream)
    
    # Look for first external offset (could be mz or intensity)
    readuntil(stream, "external offset")
    current_line = readline(stream)
    
    # Check which array comes first by looking at the param group reference
    mz_first = occursin("mzArray", current_line) ? 3 : 4
    
    # Find array length - look for the first external array length after coordinates
    seek(stream, pos_before)
    readuntil(stream, "IMS:1000103")
    array_len_skip = 0
    
    # Find spectrum end
    readuntil(stream, "</spectrum>")
    spectrum_end_skip = 0
    
    return [x_skip, y_skip, mz_first, array_len_skip, spectrum_end_skip]
end

function determine_parser(stream::IO, mz_is_compressed::Bool, int_is_compressed::Bool)
    start_pos = position(stream)
    spectrum_xml = ""
    
    try
        # Find the start of the first spectrum tag
        while !eof(stream)
            line = readline(stream)
            if occursin("<spectrum ", line)
                spectrum_buffer = IOBuffer()
                write(spectrum_buffer, line)
                # Read until the end of the spectrum tag
                while !eof(stream)
                    line = readline(stream)
                    write(spectrum_buffer, line)
                    if occursin("</spectrum>", line)
                        break
                    end
                end
                spectrum_xml = String(take!(spectrum_buffer))
                break # Found the first spectrum, so we can stop
            end
        end
    finally
        seek(stream, start_pos) # Always reset stream position
    end

    if isempty(spectrum_xml)
        # Fallback based on compression flags if no spectrum tag found
        return (mz_is_compressed || int_is_compressed) ? :compressed : :uncompressed
    end

    # Inspect the XML content
    has_neofx_markers = occursin("encodedLength=\"0\"", spectrum_xml) &&
                        occursin("external encoded length", spectrum_xml)
    
    has_external_data_markers = occursin("IMS:1000101", spectrum_xml) && 
                                occursin("IMS:1000102", spectrum_xml) &&
                                occursin("IMS:1000103", spectrum_xml)

    if has_neofx_markers
        return :neofx
    end
    
    if mz_is_compressed || int_is_compressed || has_external_data_markers
        return :compressed
    end

    return :uncompressed
end

function load_imzml_lazy(file_path::String; cache_size::Int=100)
    println("DEBUG: Checking for .imzML file at $file_path")
    if !isfile(file_path)
        throw(FileFormatError("Provided path is not a file: $(file_path)"))
    end

    ibd_path = replace(file_path, r"\.(imzML|mzML)"i => ".ibd")
    println("DEBUG: Checking for .ibd file at $ibd_path")
    if !isfile(ibd_path)
        throw(FileFormatError("Corresponding .ibd file not found for: $(file_path)"))
    end

    println("DEBUG: Opening file streams for .imzML and .ibd")
    stream = open(file_path, "r")
    hIbd = open(ibd_path, "r")

    try
        println("DEBUG: Configuring axes...")
        param_groups = axes_config_img(stream)
        println("DEBUG: Getting image dimensions...")
        imgDim = get_img_dimensions(stream)
        width, height, num_spectra = imgDim
        println("DEBUG: Image dimensions: $(width)x$(height), $num_spectra spectra.")

        # Extract default formats from the parsed param_groups
        mz_group = nothing
        int_group = nothing

        for group in values(param_groups)
            if group.Axis == 1
                mz_group = group
            elseif group.Axis == 2
                int_group = group
            end
        end

        if mz_group === nothing || int_group === nothing
            @warn "Could not find global definitions for m/z and intensity arrays. Using hardcoded defaults (Float64)."
            default_mz_format = Float64
            default_intensity_format = Float64
            mz_is_compressed = false
            int_is_compressed = false
            global_mode = UNKNOWN
        else
            default_mz_format = mz_group.Format
            default_intensity_format = int_group.Format
            mz_is_compressed = mz_group.Packed
            int_is_compressed = int_group.Packed
            global_mode = mz_group.Mode != UNKNOWN ? mz_group.Mode : int_group.Mode
        end

        println("DEBUG: m/z format: $default_mz_format, Intensity format: $default_intensity_format")
        println("DEBUG: m/z compressed: $mz_is_compressed, Intensity compressed: $int_is_compressed")
        println("DEBUG: Global mode: $global_mode")

        # --- Parser Selection ---
        parser_type = determine_parser(stream, mz_is_compressed, int_is_compressed)
        println("DEBUG: Selected parser: $parser_type")

        local spectra_metadata
        if parser_type == :neofx
            println("DEBUG: Using neofx parser.")
            spectra_metadata = parse_neofx(stream, hIbd, param_groups, width, height, num_spectra, 
                                           default_mz_format, default_intensity_format, 
                                           mz_is_compressed, int_is_compressed, global_mode)
        #=
        elseif parser_type == :compressed
            println("DEBUG: Using compressed parser.")
            spectra_metadata = parse_compressed(stream, hIbd, param_groups, width, height, num_spectra,
                                           default_mz_format, default_intensity_format,
                                           mz_is_compressed, int_is_compressed, global_mode)
        else # :uncompressed
            println("DEBUG: Using uncompressed parser.")
            spectra_metadata = parse_uncompressed(stream, hIbd, param_groups, width, height, num_spectra,
                                           default_mz_format, default_intensity_format,
                                           mz_is_compressed, int_is_compressed, global_mode)
        end
        =#
        else
            println("DEBUG: Using compressed parser.")
            spectra_metadata = parse_compressed(stream, hIbd, param_groups, width, height, num_spectra,
                                           default_mz_format, default_intensity_format,
                                           mz_is_compressed, int_is_compressed, global_mode)
        end

        println("DEBUG: Metadata parsing complete.")
        
        # Build coordinate map for imzML files
        println("DEBUG: Building coordinate map...")
        coordinate_map = zeros(Int, width, height)
        for (idx, meta) in enumerate(spectra_metadata)
            if idx == 1
                println("DIAGNOSTIC_WRITE: For index 1, attempting to write to coordinate_map[$(meta.x), $(meta.y)]")
            end
            if 1 <= meta.x <= width && 1 <= meta.y <= height
                coordinate_map[meta.x, meta.y] = idx
            end
        end
        println("DEBUG: Coordinate map built.")

        close(stream)

        source = ImzMLSource(hIbd, default_mz_format, default_intensity_format)
        println("DEBUG: Creating MSIData object.")
        return MSIData(source, spectra_metadata, (width, height), coordinate_map, cache_size)

    catch e
        close(stream)
        close(hIbd)
        rethrow(e)
    end
end

function parse_uncompressed(stream::IO, hIbd::IO, param_groups::Dict{String, SpecDim}, 
                           width::Int32, height::Int32, num_spectra::Int32,
                           mz_format::DataType, intensity_format::DataType, 
                           mz_is_compressed::Bool, int_is_compressed::Bool, global_mode::SpectrumMode)
    
    println("DEBUG: Learning file structure from first spectrum...")
    start_of_spectra_xml = position(stream)
    
    # Skip UUID at beginning of IBD file (from converter)
    seek(hIbd, 16)
    current_ibd_offset = position(hIbd)
    
    attr = get_spectrum_attributes_v2(stream, hIbd)
    seek(stream, start_of_spectra_xml)
    println("DEBUG: Initial IBD offset (after UUID): $current_ibd_offset")

    spectra_metadata = Vector{SpectrumMetadata}(undef, num_spectra)
    mz_is_first = attr[3] == 3

    for k in 1:num_spectra
        spectrum_start_pos = position(stream)
        
        # Find X coordinate
        line = readuntil(stream, "IMS:1000050")
        if eof(stream)
            error("Unexpected EOF while looking for X coordinate")
        end
        val_tag_x = find_tag(stream, r"value=\"(\d+)\"")
        x = parse(Int32, val_tag_x.captures[1])

        # Find Y coordinate  
        line = readuntil(stream, "IMS:1000051")
        val_tag_y = find_tag(stream, r"value=\"(\d+)\"")
        y = parse(Int32, val_tag_y.captures[1])

        # Find array length - look for external array length
        line = readuntil(stream, "IMS:1000103")
        val_tag_len = find_tag(stream, r"value=\"(\d+)\"")
        nPoints = parse(Int32, val_tag_len.captures[1])

        # CRITICAL FIX: Use the actual formats from the XML, not hardcoded values
        # The converter writes Float64 for mz, Float32 for intensity
        mz_len_bytes = nPoints * sizeof(Float64)  # Converter uses Float64 for mz
        int_len_bytes = nPoints * sizeof(Float32) # Converter uses Float32 for intensity

        local mz_offset, int_offset
        if mz_is_first
            mz_offset = current_ibd_offset
            int_offset = mz_offset + mz_len_bytes
        else
            int_offset = current_ibd_offset
            mz_offset = int_offset + int_len_bytes
        end

        # Mode detection
        current_pos = position(stream)
        seek(stream, spectrum_start_pos)
        spectrum_buffer = IOBuffer()
        line = ""
        while !eof(stream)
            line = readline(stream)
            write(spectrum_buffer, line)
            if occursin("</spectrum>", line)
                break
            end
        end
        spectrum_xml = String(take!(spectrum_buffer))
        
        spectrum_mode = global_mode
        if occursin("MS:1000127", spectrum_xml) 
            spectrum_mode = CENTROID
        elseif occursin("MS:1000128", spectrum_xml)
            spectrum_mode = PROFILE
        end
        seek(stream, current_pos)

        # CRITICAL FIX: Use correct data types that match the converter output
        mz_asset = SpectrumAsset(Float64, mz_is_compressed, mz_offset, nPoints, :mz)
        int_asset = SpectrumAsset(Float32, int_is_compressed, int_offset, nPoints, :intensity)
        
        spectra_metadata[k] = SpectrumMetadata(x, y, "", spectrum_mode, mz_asset, int_asset)
        
        current_ibd_offset += mz_len_bytes + int_len_bytes
        
        # Skip to next spectrum
        line = readuntil(stream, "</spectrum>")
    end

    return spectra_metadata
end

function parse_compressed(stream::IO, hIbd::IO, param_groups::Dict{String, SpecDim},
                         width::Int32, height::Int32, num_spectra::Int32,
                         default_mz_format::DataType, default_intensity_format::DataType, 
                         mz_is_compressed::Bool, int_is_compressed::Bool, global_mode::SpectrumMode)
    spectra_metadata = Vector{SpectrumMetadata}(undef, num_spectra)
    
    # Use concrete type for array_data
    array_data_type = @NamedTuple{is_mz::Bool, array_length::Int32, encoded_length::Int64, offset::Int64}

    for k in 1:num_spectra
        # Read the full spectrum XML block
        spectrum_buffer = IOBuffer()
        line = ""
        while !eof(stream)
            line = readline(stream)
            if occursin("<spectrum ", line)
                write(spectrum_buffer, line)
                break
            end
        end
        while !eof(stream)
            line = readline(stream)
            write(spectrum_buffer, line)
            if occursin("</spectrum>", line)
                break
            end
        end
        spectrum_xml = String(take!(spectrum_buffer))

        # Parse coordinates
        x_match = match(r"IMS:1000050.*?value=\"(\d+)\"", spectrum_xml)
        y_match = match(r"IMS:1000051.*?value=\"(\d+)\"", spectrum_xml)
        x = x_match !== nothing ? parse(Int32, x_match.captures[1]) : Int32(0)
        y = y_match !== nothing ? parse(Int32, y_match.captures[1]) : Int32(0)

        # Parse mode - use SpectrumMode type directly
        spectrum_mode = global_mode
        if occursin("MS:1000127", spectrum_xml) 
            spectrum_mode = CENTROID
        elseif occursin("MS:1000128", spectrum_xml)
            spectrum_mode = PROFILE
        end

        # Parse binary data arrays with concrete typing
        array_data = array_data_type[]
        
        # Find all binaryDataArray blocks
        array_matches = eachmatch(r"<binaryDataArray.*?<\/binaryDataArray>"s, spectrum_xml)
        for array_match in array_matches
            array_xml = array_match.match
            
            # Determine if this is m/z or intensity array
            is_mz = occursin("MS:1000514", array_xml) || occursin("mzArray", array_xml)
            
            # Parse external data parameters
            # Get array_length (nPoints)
            array_len_cv_match = match(r"IMS:1000103.*?value=\"(\d+)\"", array_xml)
            array_length = Int32(0)
            if array_len_cv_match !== nothing
                array_length = parse(Int32, array_len_cv_match.captures[1])
            end
            if array_length == 0
                nPoints_match = match(r"defaultArrayLength=\"(\d+)\"", spectrum_xml)
                if nPoints_match !== nothing
                    array_length = parse(Int32, nPoints_match.captures[1])
                end
            end

            # Get encoded_length
            encoded_len_cv_match = match(r"IMS:1000104.*?value=\"(\d+)\"", array_xml)
            encoded_length = Int64(0)
            if encoded_len_cv_match !== nothing
                encoded_length = parse(Int64, encoded_len_cv_match.captures[1])
            else
                encoded_len_attr_match = match(r"encodedLength=\"(\d+)\"", array_xml)
                if encoded_len_attr_match !== nothing
                    encoded_length = parse(Int64, encoded_len_attr_match.captures[1])
                end
            end

            # Get offset
            offset_match = match(r"IMS:1000102.*?value=\"(\d+)\"", array_xml)
            offset = Int64(0)
            if offset_match !== nothing
                offset = parse(Int64, offset_match.captures[1])
            end
            
            if array_length > 0 && offset > 0
                push!(array_data, (is_mz=is_mz, array_length=array_length, 
                                 encoded_length=encoded_length, offset=offset))
            end
        end

        # Separate m/z and intensity arrays
        mz_data = filter(d -> d.is_mz, array_data)
        int_data = filter(d -> !d.is_mz, array_data)
        
        # FIX: Handle empty spectra gracefully instead of throwing errors
        current_ibd_pos = position(hIbd)
        
        if length(mz_data) != 1 || length(int_data) != 1
            # Create placeholder metadata for empty/invalid spectrum
            println("DEBUG: Spectrum $k is empty or invalid - creating placeholder metadata")
            mz_asset = SpectrumAsset(default_mz_format, mz_is_compressed, current_ibd_pos, 0, :mz)
            int_asset = SpectrumAsset(default_intensity_format, int_is_compressed, current_ibd_pos, 0, :intensity)
        else
            mz_info = mz_data[1]
            int_info = int_data[1]

            # DEBUG: Print first spectrum details
            if k == 1
                println("DEBUG First spectrum parsed:")
                println("  Coordinates: x=$x, y=$y")
                println("  Mode: $spectrum_mode")
                println("  m/z array: array_length=$(mz_info.array_length), encoded_length=$(mz_info.encoded_length), offset=$(mz_info.offset)")
                println("  intensity array: array_length=$(int_info.array_length), encoded_length=$(int_info.encoded_length), offset=$(int_info.offset)")
                println("  Expected m/z bytes: $(mz_info.array_length * sizeof(default_mz_format))")
                println("  Expected intensity bytes: $(int_info.array_length * sizeof(default_intensity_format))")
            end

            # Create SpectrumAsset objects
            mz_asset = SpectrumAsset(default_mz_format, mz_is_compressed, mz_info.offset, 
                                    mz_is_compressed ? mz_info.encoded_length : mz_info.array_length, :mz)
            int_asset = SpectrumAsset(default_intensity_format, int_is_compressed, int_info.offset,
                                     int_is_compressed ? int_info.encoded_length : int_info.array_length, :intensity)
        end
        
        spectra_metadata[k] = SpectrumMetadata(x, y, "", spectrum_mode, mz_asset, int_asset)
    end

    return spectra_metadata
end

function parse_neofx(stream::IO, hIbd::IO, param_groups::Dict{String, SpecDim},
                              width::Int32, height::Int32, num_spectra::Int32,
                              default_mz_format::DataType, default_intensity_format::DataType, 
                              mz_is_compressed::Bool, int_is_compressed::Bool, global_mode::SpectrumMode)
    spectra_metadata = Vector{SpectrumMetadata}(undef, num_spectra)

    array_data_type = @NamedTuple{is_mz::Bool, array_length::Int32, encoded_length::Int64, offset::Int64}

    for k in 1:num_spectra
        # Read the full spectrum XML block
        spectrum_buffer = IOBuffer()
        line = ""
        while !eof(stream)
            line = readline(stream)
            if occursin("<spectrum ", line)
                write(spectrum_buffer, line)
                break
            end
        end
        while !eof(stream)
            line = readline(stream)
            write(spectrum_buffer, line)
            if occursin("</spectrum>", line)
                break
            end
        end
        spectrum_xml = String(take!(spectrum_buffer))

        # Parse coordinates
        x_match = match(r"IMS:1000050.*?value=\"(\d+)\"", spectrum_xml)
        y_match = match(r"IMS:1000051.*?value=\"(\d+)\"", spectrum_xml)
        x = x_match !== nothing ? parse(Int32, x_match.captures[1]) : Int32(0)
        y = y_match !== nothing ? parse(Int32, y_match.captures[1]) : Int32(0)

        # Parse mode
        spectrum_mode = global_mode
        if occursin("MS:1000127", spectrum_xml) 
            spectrum_mode = CENTROID
        elseif occursin("MS:1000128", spectrum_xml)
            spectrum_mode = PROFILE
        end

        # Parse binary data arrays
        array_data = array_data_type[]
        
        # Find all binaryDataArray blocks
        array_matches = eachmatch(r"<binaryDataArray.*?<\/binaryDataArray>"s, spectrum_xml)
        for array_match in array_matches
            array_xml = array_match.match
            
            # Determine if this is m/z or intensity array
            is_mz = occursin("MS:1000514", array_xml) || occursin("mzArray", array_xml)
            
            # Parse external data parameters
            # Get array_length (nPoints)
            array_len_cv_match = match(r"IMS:1000103.*?value=\"(\d+)\"", array_xml)
            array_length = Int32(0)
            if array_len_cv_match !== nothing
                array_length = parse(Int32, array_len_cv_match.captures[1])
            end
            if array_length == 0
                nPoints_match = match(r"defaultArrayLength=\"(\d+)\"", spectrum_xml)
                if nPoints_match !== nothing
                    array_length = parse(Int32, nPoints_match.captures[1])
                end
            end

            # Get encoded_length
            encoded_len_cv_match = match(r"IMS:1000104.*?value=\"(\d+)\"", array_xml)
            encoded_length = Int64(0)
            if encoded_len_cv_match !== nothing
                encoded_length = parse(Int64, encoded_len_cv_match.captures[1])
            else
                encoded_len_attr_match = match(r"encodedLength=\"(\d+)\"", array_xml)
                if encoded_len_attr_match !== nothing
                    encoded_length = parse(Int64, encoded_len_attr_match.captures[1])
                end
            end

            # Get offset
            offset_match = match(r"IMS:1000102.*?value=\"(\d+)\"", array_xml)
            offset = Int64(0)
            if offset_match !== nothing
                offset = parse(Int64, offset_match.captures[1])
            end
            
            if array_length > 0 && offset > 0
                push!(array_data, (is_mz=is_mz, array_length=array_length, 
                                 encoded_length=encoded_length, offset=offset))
            end
        end

        # Separate m/z and intensity arrays
        mz_data = filter(d -> d.is_mz, array_data)
        int_data = filter(d -> !d.is_mz, array_data)
        
        # FIX: Handle empty spectra gracefully instead of throwing errors
        current_ibd_pos = position(hIbd)
        
        if length(mz_data) != 1 || length(int_data) != 1
            # Create placeholder metadata for empty/invalid spectrum
            println("DEBUG: Spectrum $k is empty or invalid - creating placeholder metadata")
            mz_asset = SpectrumAsset(default_mz_format, mz_is_compressed, current_ibd_pos, 0, :mz)
            int_asset = SpectrumAsset(default_intensity_format, int_is_compressed, current_ibd_pos, 0, :intensity)
        else
            mz_info = mz_data[1]
            int_info = int_data[1]

            # DEBUG: Print first spectrum details
            if k == 1
                println("DEBUG First spectrum parsed:")
                println("  Coordinates: x=$x, y=$y")
                println("  Mode: $spectrum_mode")
                println("  m/z array: array_length=$(mz_info.array_length), encoded_length=$(mz_info.encoded_length), offset=$(mz_info.offset)")
                println("  intensity array: array_length=$(int_info.array_length), encoded_length=$(int_info.encoded_length), offset=$(int_info.offset)")
                println("  Expected m/z bytes: $(mz_info.array_length * sizeof(default_mz_format))")
                println("  Expected intensity bytes: $(int_info.array_length * sizeof(default_intensity_format))")
            end

            # Create SpectrumAsset objects
            mz_asset = SpectrumAsset(default_mz_format, mz_is_compressed, mz_info.offset, 
                                    mz_is_compressed ? mz_info.encoded_length : mz_info.array_length, :mz)
            int_asset = SpectrumAsset(default_intensity_format, int_is_compressed, int_info.offset,
                                     int_is_compressed ? int_info.encoded_length : int_info.array_length, :intensity)
        end
        
        spectra_metadata[k] = SpectrumMetadata(x, y, "", spectrum_mode, mz_asset, int_asset)
    end

    return spectra_metadata
end

# --- End of content from imzML.jl ---

# =============================================================================
#
# Image Slice Extraction
#
# =============================================================================

"""
    find_mass(mz_array, intensity_array, target_mass, tolerance)

Finds the intensity of the most intense peak within a mass tolerance window.
This optimized version uses binary search for efficiency.

# Returns
- The intensity (`Float64`) of the peak if found, otherwise `0.0`.
"""
function find_mass(mz_array::AbstractVector{<:Real}, intensity_array::AbstractVector{<:Real}, 
                   target_mass::Real, tolerance::Real)
    # Fast-path rejection: if the array is empty or the target is out of range
    if isempty(mz_array) || target_mass + tolerance < first(mz_array) || target_mass - tolerance > last(mz_array)
        return 0.0
    end

    lower_bound = target_mass - tolerance
    upper_bound = target_mass + tolerance

    # Use binary search to find the start and end of the m/z window
    start_idx = searchsortedfirst(mz_array, lower_bound)
    end_idx = searchsortedlast(mz_array, upper_bound)

    # If the window is empty, return 0.0
    if start_idx > end_idx
        return 0.0
    end

    # Find the maximum intensity within the identified window, optimized with @inbounds and @simd
    max_intensity = intensity_array[start_idx]
    @inbounds @simd for i in (start_idx + 1):end_idx
        max_intensity = max(max_intensity, intensity_array[i])
    end
    
    return max_intensity
end

"""
    load_slices(folder, masses, tolerance)

Loads image slices for multiple masses from all `.imzML` files in a directory.
This function is now refactored to use the new MSIData architecture and its
caching capabilities.
"""
function load_slices(folder::String, masses::AbstractVector{<:Real}, tolerance::Real)
    files = filter(f -> endswith(f, ".imzML"), readdir(folder, join=true))
    if isempty(files)
        @warn "No .imzML files found in the specified directory: $folder"
        return (Array{Matrix{Float64}}(undef, 0, 0), String[])
    end
    n_files = length(files)
    n_slices = length(masses)

    # FIXED: Use concrete type instead of Any
    img_list = Array{Matrix{Float64}}(undef, n_files, n_slices)
    names = String[]

    for (i, file) in enumerate(files)
        name = basename(file)
        push!(names, name)
        @info "Processing $(i)/$(n_files): $(name)"
        
        # Load data using the new lazy loader, returning an MSIData object
        msi_data = @time load_imzml_lazy(file)

        # Create empty images for all slices for the current file
        width, height = msi_data.image_dims
        current_file_slices = [zeros(Float64, width, height) for _ in 1:n_slices]

        # Use the high-performance iterator to process all spectra
        _iterate_spectra_fast(msi_data) do spec_idx, mz_array, intensity_array
            meta = msi_data.spectra_metadata[spec_idx]
            
            # Now, check for all masses of interest in this single spectrum
            for (j, mass) in enumerate(masses)
                intensity = find_mass(mz_array, intensity_array, mass, tolerance)
                if intensity > 0.0
                    if 1 <= meta.x <= width && 1 <= meta.y <= height
                        current_file_slices[j][meta.y, meta.x] = intensity
                    end
                end
            end
        end # end of fast iterator

        # Assign the generated images to the main list
        for j in 1:n_slices
            img_list[i, j] = current_file_slices[j]
        end

    end # end of files loop
    
    return (img_list, names)
end


"""
    get_mz_slice(data::MSIData, mass::Real, tolerance::Real)

Extracts an image slice for a given m/z value without plotting.
This is a performant function that iterates through spectra once.

# Returns
- A `Matrix{Float64}` representing the intensity slice.
"""
function get_mz_slice(data::MSIData, mass::Real, tolerance::Real; mask_path::Union{String, Nothing}=nothing)
    width, height = data.image_dims
    slice_matrix = zeros(Float64, height, width)

    local masked_indices::Union{Set{Int}, Nothing} = nothing
    if mask_path !== nothing
        mask_matrix = load_and_prepare_mask(mask_path, (width, height))
        masked_indices = get_masked_spectrum_indices(data, mask_matrix)
        println("Applying mask from: $(mask_path), found $(length(masked_indices)) spectra in ROI.")
    end

    # INTELLIGENT LOADING: Ensure analytics are computed for filtering.
    if data.spectrum_stats_df === nothing || !hasproperty(data.spectrum_stats_df, :MinMZ)
        println("Per-spectrum metadata not found. Running one-time analytics computation...")
        precompute_analytics(data)
    end

    println("Using high-performance sequential iterator...")
    target_min = mass - tolerance
    target_max = mass + tolerance
    stats_df = data.spectrum_stats_df

    # 1. Find all candidate spectra first for efficient filtering
    candidate_indices = Set{Int}()
    indices_to_check = masked_indices === nothing ? (1:length(data.spectra_metadata)) : masked_indices
    
    for i in indices_to_check
        spec_min_mz = stats_df.MinMZ[i]
        spec_max_mz = stats_df.MaxMZ[i]
        if target_max >= spec_min_mz && target_min <= spec_max_mz
            push!(candidate_indices, i)
        end
    end

    println("Found $(length(candidate_indices)) candidate spectra (filtered from $(length(indices_to_check)) initial spectra)")

    # 2. Iterate using the optimized, low-allocation iterator
    results_count = 0
    _iterate_spectra_fast(data, collect(candidate_indices)) do idx, mz_array, intensity_array
        meta = data.spectra_metadata[idx]
        intensity = find_mass(mz_array, intensity_array, mass, tolerance)
        if intensity > 0.0
            if 1 <= meta.x <= width && 1 <= meta.y <= height
                slice_matrix[meta.y, meta.x] = intensity
                results_count += 1
            end
        end
    end
    
    println("Populated $results_count pixels with intensity data")
    replace!(slice_matrix, NaN => 0.0)
    return slice_matrix
end


"""
    get_multiple_mz_slices(data::MSIData, masses::AbstractVector{<:Real}, tolerance::Real)

Extracts multiple image slices for a given list of m/z values in a single pass.
This is a highly performant function that iterates through the full dataset only once.

# Returns
- A `Dict{Real, Matrix{Float64}}` mapping each mass to its intensity slice matrix.
"""
function get_multiple_mz_slices(data::MSIData, masses::AbstractVector{<:Real}, tolerance::Real; mask_path::Union{String, Nothing}=nothing)
    width, height = data.image_dims
    
    # Sort masses to improve cache locality during search
    sorted_masses = sort(masses)

    # 1. Initialize a dictionary to hold the output slice matrices
    slice_dict = Dict{Real, Matrix{Float64}}()
    for mass in sorted_masses
        slice_dict[mass] = zeros(Float64, height, width)
    end

    local masked_indices::Union{Set{Int}, Nothing} = nothing
    if mask_path !== nothing
        mask_matrix = load_and_prepare_mask(mask_path, (width, height))
        masked_indices = get_masked_spectrum_indices(data, mask_matrix)
        println("Applying mask from: $(mask_path), found $(length(masked_indices)) spectra in ROI.")
    end

    # 2. Ensure analytics are computed for filtering.
    if data.spectrum_stats_df === nothing || !hasproperty(data.spectrum_stats_df, :MinMZ)
        println("Per-spectrum metadata not found. Running one-time analytics computation...")
        precompute_analytics(data)
    end

    println("Filtering candidate spectra for $(length(masses)) m/z values...")
    stats_df = data.spectrum_stats_df
    candidate_indices = Set{Int}()
    indices_to_check = masked_indices === nothing ? (1:length(data.spectra_metadata)) : masked_indices

    # 3. Find all spectra that could contain *any* of the requested masses.
    for mass in sorted_masses
        target_min = mass - tolerance
        target_max = mass + tolerance
        for i in indices_to_check
            # If already a candidate, no need to check again
            if i in candidate_indices
                continue
            end
            spec_min_mz = stats_df.MinMZ[i]
            spec_max_mz = stats_df.MaxMZ[i]
            if target_max >= spec_min_mz && target_min <= spec_max_mz
                push!(candidate_indices, i)
            end
        end
    end

    println("Found $(length(candidate_indices)) total candidate spectra.")

    # 4. Iterate through the data a single time using the optimized iterator.
    _iterate_spectra_fast(data, collect(candidate_indices)) do idx, mz_array, intensity_array
        meta = data.spectra_metadata[idx]
        # For this single spectrum, check all masses of interest
        for mass in sorted_masses
            # Check if this spectrum's range actually covers the current mass
            # This is a finer-grained check than the initial filtering
            if !isempty(mz_array) && (mass + tolerance) >= first(mz_array) && (mass - tolerance) <= last(mz_array)
                intensity = find_mass(mz_array, intensity_array, mass, tolerance)
                if intensity > 0.0
                    if 1 <= meta.x <= width && 1 <= meta.y <= height
                        slice_dict[mass][meta.y, meta.x] = intensity
                    end
                end
            end
        end
    end
    
    # 5. Clean up and return
    for mass in sorted_masses
        replace!(slice_dict[mass], NaN => 0.0)
    end
    
    println("Finished generating $(length(masses)) slices in a single pass.")
    return slice_dict
end


"""
    plot_slice(msi_data::MSIData, mass::Float64, tolerance::Float64, output_dir::String; stage_name="slice", bins=256)

Generates and saves a plot of a single image slice for a given m/z value.
This function closely imitates the logic of the original `GetSlice` but uses
the modern `MSIData` access patterns and robust peak finding.
"""
function plot_slice(msi_data::MSIData, mass::Real, tolerance::Real, output_dir::String; 
                    stage_name="slice_mz_$(mass)", bins=256, mask_path::Union{String, Nothing}=nothing)
    
    # 1. Generate the slice, applying the mask if provided.
    println("Generating slice for m/z $mass...")
    slice_matrix = get_mz_slice(msi_data, mass, tolerance, mask_path=mask_path)
    println("Slice generation complete.")

    # 2. Plot the resulting slice matrix using CairoMakie
    println("Plotting slice...")
    
    fig = Figure(size = (600, 500))
    ax = CairoMakie.Axis(fig[1, 1],
        aspect=DataAspect(),
        title=@sprintf("Slice for m/z: %.2f", mass),
        yreversed=true
    )
    hidedecorations!(ax)

    # Use mass-specific bounds for colorrange, ensuring a valid range
    min_val, max_val = extrema(slice_matrix)
    if min_val == max_val
        max_val = min_val + 1.0 # Ensure the color range has a non-zero width
    end
    hm = heatmap!(ax, slice_matrix,
        colormap=cgrad(ColorSchemes.viridis, bins),
        colorrange=(min_val, max_val)
    )
    
    Colorbar(fig[1, 2], hm, label="Intensity")
    colgap!(fig.layout, 5)

    # 3. Save the plot
    mkpath(output_dir)
    filename = "$(stage_name).png"
    save_path = joinpath(output_dir, filename)
    save(save_path, fig)
    @info "Saved slice plot to $save_path"
    
    return fig
end

# ============================================================================
#
#
# Image Processing and Normalization
#
# ============================================================================

"""
    get_outlier_thres(img, prob=0.98)

Calculates dynamic intensity range bounds for an image based on a cumulative
probability histogram. This function replicates an R algorithm to find an
intensity threshold that corresponds to a given cumulative probability, which
is used to exclude outliers before normalization.

# Arguments
- `img`: The input image matrix.
- `prob`: The cumulative probability threshold (default: 0.98) for outlier detection.

# Returns
- A tuple `(low, high)` representing the calculated lower and upper intensity bounds.
"""
function get_outlier_thres(img::AbstractVector{<:Real}, prob::Real=0.98)
    # DO NOT filter zeros. Use all pixel values like R does.
    int_values = img
    low = minimum(int_values)
    upp = maximum(int_values)

    # Create histogram bins
    if upp - low + 1 >= 100
        bins = 100
        brk = range(low, stop=upp + (upp - low)/(bins - 1), length=bins + 1)
    else
        brk = collect(floor(low):(ceil(upp) + 1))
    end

    # Compute histogram with RIGHT-CLOSED = FALSE to mimic R's right=FALSE
    # This makes intervals [a, b) instead of the default [a, b]
    h = fit(Histogram, int_values, brk, closed=:left) # Key change: closed=:left

    # Check if histogram is valid
    if isempty(h.weights) || sum(h.weights) == 0
        return (low, upp)
    end

    # Replicate R's algorithm exactly
    cum_counts = cumsum(h.weights) / sum(h.weights)
    # Calculate the difference from the target probability
    top = prob .- cum_counts
    delta = abs.(top)
    # Find the index where the difference is minimized
    min_delta_index = findfirst(x -> x == minimum(delta), delta)
    
    # R adds 1 to this index: index <- 1 + which(...)[1]
    target_bin_index = min_delta_index + 1
    # Get the upper edge of the target bin
    target_bin_upper_edge = h.edges[1][target_bin_index]

    # Find the maximum data value that is strictly less than this edge
    # This mimics: max( intMap[ intMap < h$breaks[index] ] )
    values_below_edge = filter(x -> x < target_bin_upper_edge, int_values)
    actual_threshold = isempty(values_below_edge) ? low : maximum(values_below_edge)

    return (low, actual_threshold)
end

function get_outlier_thres(img::AbstractMatrix{<:Real}, prob::Real=0.98)
    return get_outlier_thres(vec(img), prob)
end

"""
    set_pixel_depth(img, bounds, depth)

Quantizes the intensity values of an image into a specified number of bins
(`depth`) within a given intensity range (`bounds`). Pixels outside the bounds
are clipped.

# Arguments
- `img`: The input image matrix.
- `bounds`: A tuple `(min, max)` specifying the intensity range for quantization.
- `depth`: The number of quantization levels (bins).

# Returns
- A `Matrix{UInt8}` with pixel values quantized to the specified depth.
"""
function set_pixel_depth(img::AbstractMatrix{<:Real}, bounds::Tuple{<:Real, <:Real}, depth::Integer)
    min_val, max_val = bounds
    bins = depth - 1
    
    if min_val >= max_val
        return zeros(UInt8, size(img))
    end
    
    # Create intensity bins - FIXED: use binary search instead of linear scan
    range_vals = range(min_val, stop=max_val, length=depth)[2:depth]
    
    # Pre-allocate result for better performance
    result = similar(img, UInt8)
    
    # Use binary search for much faster bin assignment
    @inbounds for i in eachindex(img)
        val = img[i]
        if val <= min_val
            result[i] = 0
        elseif val >= max_val
            result[i] = bins
        else
            # Binary search is much faster than linear scan
            bin_idx = searchsortedfirst(range_vals, val)
            result[i] = bin_idx
        end
    end
    
    return result
end

"""
    TrIQ(pixMap, depth, prob=0.98)

Applies TrIQ (Treshold Intensity Quantization) normalization to an image.
This function first computes dynamic intensity range bounds by identifying outliers
based on a cumulative probability, then sets the pixel depth (quantizes intensities)
within these bounds.

# Arguments
- `pixMap`: The input image matrix (e.g., a slice from `GetMzSliceJl`).
- `depth`: The number of grey levels (bins) to quantize the intensities into.
- `prob`: The target cumulative probability (e.g., 0.98) used to determine outlier thresholds.

# Returns
- A new image matrix with intensities quantized to the specified depth within the TrIQ bounds.
"""
function TrIQ(pixMap::AbstractMatrix{<:Real}, depth::Integer, prob::Real=0.98; mask_matrix::Union{BitMatrix, Nothing}=nothing)
    local values_for_thres
    if mask_matrix !== nothing
        values_for_thres = pixMap[mask_matrix]
        filter!(x -> x != 0, values_for_thres)
    else
        values_for_thres = pixMap
    end

    if isempty(values_for_thres) || all(iszero, values_for_thres)
        bounds = (0.0, 0.0)
        quantized = zeros(UInt8, size(pixMap))
    else
        bounds = get_outlier_thres(values_for_thres, prob)
        quantized = set_pixel_depth(pixMap, bounds, depth)
    end
    
    return quantized, bounds  # Return both!
end

"""
    norm_slices_hist(slices, bins; prob=0.98)

Normalizes a set of image slices based on a shared histogram range.
"""
function norm_slices_hist(slices, bins; prob=0.98)
    n_files, n_masses = size(slices)
    norm_img = similar(slices)
    mass_bounds = []  # This will store bounds for EACH mass
    
    # Calculate bounds for each mass across all files
    for mass_idx in 1:n_masses
        # Get all slices for this specific mass across all files
        mass_slices = [slices[i, mass_idx] for i in 1:n_files]
        all_vals = reduce(vcat, [vec(s) for s in mass_slices])
        
        # Calculate global bounds for this specific mass
        mass_global_bounds = get_outlier_thres(all_vals, prob)
        push!(mass_bounds, mass_global_bounds)
        
        # Normalize each file's slice for this mass using its specific bounds
        for file_idx in 1:n_files
            norm_img[file_idx, mass_idx] = set_pixel_depth(slices[file_idx, mass_idx], mass_global_bounds, bins)
        end
    end
    
    return (norm_img=norm_img, bounds=mass_bounds)  # bounds is now a vector, one per mass
end

"""
    quantize_intensity(slice::AbstractMatrix{<:Real}, levels::Integer=256)

Linearly scales the intensity values in a slice to a specified number of levels.
The output is an array of `UInt8` values.
This is a modernized version of `IntQuantCl`.
"""
function quantize_intensity(slice::AbstractMatrix{<:Real}, levels::Integer=256; mask_matrix::Union{BitMatrix, Nothing}=nothing)
    local max_val
    if mask_matrix !== nothing
        masked_slice = slice[mask_matrix]
        # Filter out zeros that might have been introduced by masking, if any
        filter!(x -> x != 0, masked_slice)
        max_val = isempty(masked_slice) ? 0.0 : maximum(masked_slice)
    else
        max_val = maximum(slice)
    end

    if max_val <= 0
        bounds = (0.0, 0.0)
        quantized = zeros(UInt8, size(slice))
    else
        bounds = (0.0, max_val)
        quantized = similar(slice, UInt8)
        scale = (levels - 1) / max_val
        @inbounds for i in eachindex(slice)
            quantized[i] = round(UInt8, clamp(slice[i] * scale, 0, levels - 1))
        end
    end
    
    return quantized, bounds  # Return both!
end

"""
    median_filter(img)

Applies a 3x3 median filter to the input image. This is a simple noise
reduction technique that replaces each pixel's value with the median value of
its 3x3 neighborhood.

# Arguments
- `img`: The input image matrix.

# Returns
- A new matrix containing the filtered image.
"""
function median_filter(img)
    # 3x3 median filter implementation
    return mapwindow(median, img, (3, 3))
end



"""
    downsample_spectrum(mz, intensity, n_points=2000)

Reduces the number of points in a spectrum for faster plotting, while preserving peaks.
It divides the m/z range into `n_points` bins and keeps only the most intense point from each bin.
"""
function downsample_spectrum(mz::AbstractVector, intensity::AbstractVector, n_points::Integer=2000)
    if isempty(mz) || length(mz) <= n_points
        return mz, intensity
    end

    mz_min, mz_max = extrema(mz)
    bin_width = (mz_max - mz_min) / n_points
    
    # We use a vector of tuples to store the max intensity and its corresponding mz for each bin
    # (max_intensity, mz_value)
    bins = fill((0.0, 0.0), n_points)
    
    for i in eachindex(mz)
        # Determine the bin for the current point
        # Bin indices are 1-based
        bin_index = min(n_points, floor(Int, (mz[i] - mz_min) / bin_width) + 1)
        
        # If the current point's intensity is higher than what's in the bin, replace it
        if intensity[i] > bins[bin_index][1]
            bins[bin_index] = (intensity[i], mz[i])
        end
    end
    
    # Filter out empty bins and separate the mz and intensity values
    final_mz = [b[2] for b in bins if b[1] > 0.0]
    final_intensity = [b[1] for b in bins if b[1] > 0.0]
    
    return final_mz, final_intensity
end

# ============================================================================
#
#
# Analysis and Visualization
#
# ============================================================================

"""
    display_statistics(slices::AbstractArray{<:AbstractMatrix{<:Real}, 2}, 
                           names::AbstractVector{String}, masses::AbstractVector{<:Real})

Calculates and prints key statistics for each slice.
"""
function display_statistics(slices::AbstractArray{<:AbstractMatrix{<:Real}, 2}, 
                           names::AbstractVector{String}, masses::AbstractVector{<:Real})
    if isempty(slices)
        @warn "Cannot display statistics for empty slice list."
        return nothing
    end

    n_files, n_masses = size(slices)
    stats_to_calc = Dict{String, Function}("Mean" => mean)
    all_dfs = Dict{String, DataFrame}()

    for (stat_name, stat_func) in stats_to_calc
        # Pre-allocate matrix for better performance
        stat_matrix = zeros(Float64, n_files, n_masses)
        
        @inbounds for i in 1:n_files, j in 1:n_masses
            flat_slice = vec(slices[i, j])
            if !isempty(flat_slice)
                stat_matrix[i, j] = stat_func(flat_slice)
            end
        end

        # Create the DataFrame with masses as column headers
        df = DataFrame(stat_matrix, Symbol.(masses))
        # Insert the file names as the first column
        insertcols!(df, 1, :Data => names)
        mz_row = ["m/z"; masses...]
        push!(df, mz_row)

        println("\n--- Statistics: $(stat_name) ---")
        println(df)
        all_dfs[stat_name] = df
    end
    return all_dfs
end

"""
    plot_slices(slices, names, masses, output_dir; stage_name, bins=256, dpi=150, global_bounds=nothing)

Generates and saves a grid of image slice plots. Each row corresponds to a
file and each column to a mass, creating a comprehensive overview.

# Arguments
- `slices`: A 2D array of image slice matrices (`n_files` x `n_masses`).
- `names`: A vector of file names, used for titling rows.
- `masses`: A vector of m/z values, used for titling columns.
- `output_dir`: The directory where the output plots will be saved.

# Keyword Arguments
- `stage_name`: A string used to name the output files (e.g., "raw", "normalized").
- `bins`: The number of color levels in the heatmap palette.
- `dpi`: The resolution for the saved image files.
- `global_bounds`: A vector of `(min, max)` tuples, one for each mass, to ensure a consistent
  color scale across all files for a given mass. If `nothing`, bounds are calculated automatically.

# Returns
- The generated `Figure` object from Makie.
"""
function plot_slices(slices, names, masses, output_dir; stage_name, bins=256, dpi=150, global_bounds=nothing)
    n_files, n_masses = size(slices)
    
    mkpath(output_dir)

    # If global_bounds is provided, it should now be a VECTOR of bounds (one per mass)
    # If not provided, calculate per-mass bounds
    if global_bounds === nothing
        global_bounds = []
        for mass_idx in 1:n_masses
            mass_slices = [slices[i, mass_idx] for i in 1:n_files]
            all_vals = reduce(vcat, [vec(s) for s in mass_slices])
            filter!(isfinite, all_vals)
            mass_bounds = isempty(all_vals) ? (0.0, 1.0) : extrema(all_vals)
            push!(global_bounds, mass_bounds)
        end
    end

    # With Makie, we define a Figure and a layout.
    fig = Figure(size = (400 * n_masses, 330 * n_files))  # Adjusted size calculation

    # Add a title for the entire figure.
    Label(fig[0, 1:(2*n_masses)], "Image Slices - $(stage_name)", fontsize=24, font=:bold, tellwidth=false, padding=(0,0,10,0))

    for file_idx in 1:n_files
        for mass_idx in 1:n_masses
            img = slices[file_idx, mass_idx]
            mass = masses[mass_idx]
            name = names[file_idx]
            mass_global_bounds = global_bounds[mass_idx]  # Get bounds for THIS specific mass

            # --- Calculate tick properties for THIS mass ---
            min_val, max_val = mass_global_bounds
            levels = range(min_val, stop=max_val, length=bins + 1)
            level_range = levels[end] - levels[1]
            
            if level_range == 0
                levels = range(min_val - 0.1, stop=max_val + 0.1, length=bins + 1)
                level_range = 0.2
            end
            
            exponent = level_range > 0 ? floor(log10(level_range)) / 3 : 0
            scale = 10^(3 * exponent)
            scaled_levels = levels ./ scale
            
            format_num = level_range > 0 ? floor(log10(level_range)) % 3 : 0
            labels = if format_num == 0
                [@sprintf("%3.2f", lvl) for lvl in scaled_levels]
            elseif format_num == 1
                [@sprintf("%3.2f", lvl) for lvl in scaled_levels]
            else
                [@sprintf("%3.2f", lvl) for lvl in scaled_levels]
            end
            
            divisors = 2:7
            remainders = (bins - 1) .% divisors
            best_divisor = divisors[findlast(x -> x == minimum(remainders), remainders)]
            tick_indices = round.(Int, range(1, stop=bins + 1, length=best_divisor + 1))

            if !(1 in tick_indices)
                pushfirst!(tick_indices, 1)
            end
            if !((bins + 1) in tick_indices)
                push!(tick_indices, bins + 1)
            end
            unique!(sort!(tick_indices))
            
            tick_positions = levels[tick_indices]
            tick_labels = labels[tick_indices]
            # --- End of mass-specific tick calculation ---

            # Create an Axis for the heatmap
            ax = CairoMakie.Axis(fig[file_idx, 2*mass_idx-1],
                aspect=DataAspect(),
                title=@sprintf("%s\nm/z: %.2f", basename(name), mass),
                titlesize=14
            )
            hidedecorations!(ax)

            # Use mass-specific bounds for colorrange
            hm = heatmap!(ax, transpose(img),
                colormap=cgrad(ColorSchemes.viridis, bins),
                colorrange=mass_global_bounds  # ← This is mass-specific
            )

            # Add a colorbar with mass-specific scale
            cb = Colorbar(fig[file_idx, 2*mass_idx], hm,
                label=(scale == 1 ? "" : "×10^$(round(Int, 3 * exponent))"),
                labelpadding=2,
                labelsize=12,
                ticks=(tick_positions, tick_labels),
                ticklabelsize=10
            )
            colsize!(fig.layout, 2*mass_idx, 30)
        end
    end

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 10)

    # Save in multiple formats.
    formats = ["png", "pdf"]
    
    for fmt in formats
        filename = "$(stage_name)_overview.$(fmt)"
        save_path = joinpath(output_dir, filename)
        save(save_path, fig, px_per_unit = dpi / 96.0)
        @info "Saved $(fmt) overview plot to $save_path"
    end
    return fig
end

"""
    save_bitmap(name::String, pixMap::Matrix{UInt8}, colorTable::Vector{UInt32})

Saves an 8-bit indexed image as a BMP file.
This is a modernized version of `SaveBitmapCl`.
"""
function save_bitmap(name::String, pixMap::Matrix{UInt8}, colorTable::Vector{UInt32})
    # Get image dimensions
    height, width = size(pixMap)

    # Normalize pixel values to stretch contrast, as in the original SaveBitmapCl
    minVal, maxVal = extrema(pixMap)
    if maxVal > minVal
        pixMap = round.(UInt8, 255 * (pixMap .- minVal) ./ (maxVal - minVal))
    end

    # Compute row padding (each row must be a multiple of 4 bytes)
    padding = (4 - (width % 4)) % 4

    # BMP color table must have 256 entries for 8-bit images
    fullColorTable = Vector{UInt32}(undef, 256)
    if length(colorTable) <= 256
        fullColorTable[1:length(colorTable)] .= colorTable
        fullColorTable[length(colorTable)+1:end] .= 0
    else
        fullColorTable .= colorTable[1:256]
    end

    # Compute file dimensions
    offset = 14 + 40 + (256 * 4) # 14(file) + 40(info) + 1024(palette) = 1078
    imgBytes = height * (width + padding)
    fileSize = offset + imgBytes

    open(name, "w") do stream
        # === File Header (14 bytes) ===
        write(stream, UInt16(0x4D42)) # "BM"
        write(stream, UInt32(fileSize))
        write(stream, UInt16(0)) # Reserved
        write(stream, UInt16(0)) # Reserved
        write(stream, UInt32(offset))

        # === Info Header (40 bytes) ===
        write(stream, UInt32(40)) # Info header size
        write(stream, Int32(width))
        write(stream, Int32(height)) # Positive for bottom-up storage
        write(stream, UInt16(1)) # Number of color planes
        write(stream, UInt16(8)) # Bits per pixel
        write(stream, UInt32(0)) # Compression (BI_RGB)
        write(stream, UInt32(imgBytes))
        write(stream, Int32(2835)) # Pels per meter X (~72 DPI)
        write(stream, Int32(2835)) # Pels per meter Y (~72 DPI)
        write(stream, UInt32(256)) # Colors in color table
        write(stream, UInt32(0)) # Important colors (0 = all)

        # === Color Table ===
        write(stream, fullColorTable)

        # === Image Pixels (written bottom-up) ===
        row_buffer = Vector{UInt8}(undef, width + padding)
        row_buffer[width+1:end] .= 0 # pre-fill padding bytes

        for i in height:-1:1
            row_data = @view pixMap[i, :]
            row_buffer[1:width] = row_data
            write(stream, row_buffer)
        end
    end
end

"""
    generate_palette(colorscheme, n_colors=256)

Generates a BMP-compatible UInt32 color palette from a ColorScheme.
"""
function generate_palette(colorscheme, n_colors=256)
    palette = Vector{UInt32}(undef, n_colors)
    colors = get(colorscheme, range(0, 1, length=n_colors))
    for i in 1:n_colors
        c = colors[i]
        r = round(UInt8, c.r * 255)
        g = round(UInt8, c.g * 255)
        b = round(UInt8, c.b * 255)
        # BMP color table format is 0x00RRGGBB, written little-endian becomes BB GG RR 00
        palette[i] = (UInt32(r) << 16) | (UInt32(g) << 8) | UInt32(b)
    end
    return palette
end

# ============================================================================
#
#
# Additional Performance Optimizations
#
# ============================================================================

"""
    precompute_spectrum_bounds(msi_data::MSIData)

Precomputes min/max m/z values for each spectrum to enable fast filtering.
This avoids repeated computation during slice extraction.
"""
function precompute_spectrum_bounds(msi_data::MSIData)
    n_spectra = length(msi_data.spectra_metadata)
    min_mz = Vector{Float64}(undef, n_spectra)
    max_mz = Vector{Float64}(undef, n_spectra)
    
    _iterate_spectra_fast(msi_data) do idx, mz_array, intensity_array
        min_mz[idx] = first(mz_array)
        max_mz[idx] = last(mz_array)
    end
    
    return (min_mz, max_mz)
end

"""
    process_columns_first!(A::AbstractMatrix)

Process a matrix column-by-column for better cache performance.
In Julia, arrays are stored in column-major order.
"""
function process_columns_first!(A::AbstractMatrix)
    nrows, ncols = size(A)
    @inbounds for col in 1:ncols, row in 1:nrows
        # Process A[row, col] here
        # This order is ~70% faster due to cache locality
    end
    return A
end

# Viridis color palette (256 colors) - defined as constant
const ViridisPalette = generate_palette(ColorSchemes.viridis)

function generate_colorbar_image(slice_data::AbstractMatrix, color_levels::Int, output_path::String, 
                               bounds::Tuple{Float64, Float64}; 
                               use_triq::Bool=false, triq_prob::Float64=0.98, 
                               mask_path::Union{String, Nothing}=nothing)
    # Use the provided bounds instead of recalculating
    min_val, max_val = bounds
    
    # If bounds are invalid, calculate fallback
    if min_val == max_val == 0.0
        local data_for_bounds
        if mask_path !== nothing
            height, width = size(slice_data)
            mask_matrix = load_and_prepare_mask(mask_path, (width, height))
            data_for_bounds = slice_data[mask_matrix]
            filter!(x -> x > 0, data_for_bounds)
        else
            data_for_bounds = vec(slice_data)
        end
        min_val, max_val = isempty(data_for_bounds) ? (0.0, 1.0) : extrema(data_for_bounds)
    end

    # 2. Replicate the tick calculation logic from plot_slices
    bins = color_levels
    levels = range(min_val, stop=max_val, length=bins + 1)
    level_range = levels[end] - levels[1]
    
    if level_range == 0
        levels = range(min_val - 0.1, stop=max_val + 0.1, length=bins + 1)
        level_range = 0.2
    end
    
    exponent = level_range > 0 ? floor(log10(level_range)) / 3 : 0
    scale = 10^(3 * exponent)
    scaled_levels = levels ./ scale
    
    format_num = level_range > 0 ? floor(log10(level_range)) % 3 : 0
    labels = if format_num == 0
        [ @sprintf("%3.2f", lvl) for lvl in scaled_levels]
    elseif format_num == 1
        [ @sprintf("%3.2f", lvl) for lvl in scaled_levels]
    else
        [ @sprintf("%3.2f", lvl) for lvl in scaled_levels]
    end
    
    divisors = 2:7
    remainders = (bins - 1) .% divisors
    best_divisor = divisors[findlast(x -> x == minimum(remainders), remainders)]
    tick_indices = round.(Int, range(1, stop=bins + 1, length=best_divisor + 1))

    if !(1 in tick_indices)
        pushfirst!(tick_indices, 1)
    end
    if !((bins + 1) in tick_indices)
        push!(tick_indices, bins + 1)
    end
    unique!(sort!(tick_indices))
    
    tick_positions = levels[tick_indices]
    tick_labels = labels[tick_indices]

    # 3. Create and save the colorbar image
    fig = Figure(size=(150, 250))
    Colorbar(fig[1, 1], 
        colormap=cgrad(:viridis, bins, categorical=true),
        # limits=(min_val, max_val),
        limits=(levels[1], levels[end]),
        label=(scale == 1 ? "Intensity" : "Intensity ×10^$(round(Int, 3 * exponent))"),
        ticks=(tick_positions, tick_labels),
        labelsize=20,
        ticklabelsize=16
    )
    save(output_path, fig)
end
