# src/mzML.jl

# This file is responsible for parsing metadata from .mzML files.
# It has been refactored to produce a unified MSIData object.

# Constants for CV parameter accessions - defined once for performance
const MZ_AXIS_ACCESSION = "MS:1000514"
const INTENSITY_AXIS_ACCESSION = "MS:1000515"
const COMPRESSION_ACCESSION = "MS:1000574"
const NO_COMPRESSION_ACCESSION = "MS:1000576"
const CENTROID_SPECTRUM_ACCESSION = "MS:1000127"
const PROFILE_SPECTRUM_ACCESSION = "MS:1000128"
const INSTRUMENT_MODEL_ACCESSION = "MS:1000031"
const MASS_RESOLVING_POWER_ACCESSION = "MS:1001496"
const RESOLUTION_ACCESSION = "MS:1000011"
const MASS_ACCURACY_PPM_ACCESSION = "MS:1000016"
const POSITIVE_SCAN_ACCESSION = "MS:1000130"
const NEGATIVE_SCAN_ACCESSION = "MS:1000129"
const EXTERNAL_CALIBRATION_ACCESSION = "MS:1000592"
const INTERNAL_CALIBRATION_ACCESSION = "MS:1000593"
const INSTRUMENT_SPECIFIC_CALIBRATION_ACCESSION = "MS:1000747"
const LASER_WAVELENGTH_ACCESSION = "MS:1000867"
const LASER_FLUENCE_ACCESSION = "MS:1000868"
const LASER_REPETITION_RATE_ACCESSION = "MS:1000869"
const BASELINE_CORRECTION_ACCESSION = "MS:1000579"
const SMOOTHING_ACCESSION = "MS:1000580"
const DATA_TRANSFORMATION_ACCESSION = "MS:1000578"
const DEISOTOPING_ACCESSION = "MS:1000800"

# Data format accessions as constants
const DATA_FORMAT_ACCESSIONS = Dict{String, DataType}(
    "MS:1000518" => Int16,
    "MS:1000519" => Int32,
    "MS:1000520" => Float64,
    "MS:1000521" => Float32,
    "MS:1000522" => Int64,
    "MS:1000523" => Float64
)

function parse_instrument_metadata_mzml(stream::IO)
    println("DEBUG: Starting mzML instrument metadata parsing...")
    # Initialize with default values from the InstrumentMetadata constructor
    instrument_meta = InstrumentMetadata()
    
    # Create temporary variables to hold parsed values
    resolution = instrument_meta.resolution
    instrument_model = instrument_meta.instrument_model
    mass_accuracy_ppm = instrument_meta.mass_accuracy_ppm
    polarity = instrument_meta.polarity
    calibration_status = instrument_meta.calibration_status
    laser_settings = Dict{String, Any}() # Use Any for heterogeneous values
    vendor_preprocessing_steps = String[] # Initialize vendor_preprocessing_steps

    try
        # Reset stream to ensure we read from the beginning of the relevant section
        seekstart(stream)
        
        # Define state variables for parsing
        in_instrument_config_list = false
        in_data_processing_list = false

        header_read_limit = 50000 # Read up to 50KB of header, or until main data section
        bytes_read = 0
        current_pos = position(stream)

        while !eof(stream) && bytes_read < header_read_limit
            line_pos = position(stream)
            line = readline(stream)
            bytes_read += (position(stream) - line_pos)
            
            # Stop at the start of the main data section
            if occursin("<run>", line) || occursin("<spectrumList>", line)
                break
            end

            # State management for nested sections
            if occursin("<instrumentConfigurationList", line)
                in_instrument_config_list = true
            elseif occursin("</instrumentConfigurationList>", line)
                in_instrument_config_list = false
            elseif occursin("<dataProcessingList", line)
                in_data_processing_list = true
            elseif occursin("</dataProcessingList>", line)
                in_data_processing_list = false
            end
            
            if occursin("<cvParam", line)
                accession_match = match(ACCESSION_REGEX, line)
                value_match = match(VALUE_REGEX, line)
                name_match = match(NAME_REGEX, line)
                
                if accession_match !== nothing
                    acc = accession_match.captures[1]
                    val = (value_match !== nothing) ? value_match.captures[1] : ""
                    name = (name_match !== nothing) ? name_match.captures[1] : ""

                    if in_instrument_config_list # Process instrument configuration details
                        if acc == INSTRUMENT_MODEL_ACCESSION # instrument model
                            instrument_model = val
                        elseif acc == MASS_RESOLVING_POWER_ACCESSION # mass resolving power (more specific)
                            resolution = tryparse(Float64, val)
                        elseif acc == RESOLUTION_ACCESSION && resolution === nothing # resolution (less specific)
                            resolution = tryparse(Float64, val)
                        elseif acc == MASS_ACCURACY_PPM_ACCESSION # mass accuracy (ppm)
                            mass_accuracy_ppm = tryparse(Float64, val)
                        elseif acc == POSITIVE_SCAN_ACCESSION # positive scan
                            polarity = :positive
                        elseif acc == NEGATIVE_SCAN_ACCESSION # negative scan
                            polarity = :negative
                        end
                    elseif in_data_processing_list # Process vendor preprocessing steps
                        if acc == BASELINE_CORRECTION_ACCESSION # baseline correction
                            push!(vendor_preprocessing_steps, "Baseline Correction")
                        elseif acc == SMOOTHING_ACCESSION # smoothing
                            push!(vendor_preprocessing_steps, "Smoothing")
                        elseif acc == DATA_TRANSFORMATION_ACCESSION # data transformation (e.g., centroiding)
                            push!(vendor_preprocessing_steps, "Data Transformation: $(name)")
                        elseif acc == DEISOTOPING_ACCESSION # deisotoping
                            push!(vendor_preprocessing_steps, "Deisotoping")
                        end
                    else # General CV params outside specific lists
                        if acc == EXTERNAL_CALIBRATION_ACCESSION # external calibration
                            calibration_status = :external
                        elseif acc == INTERNAL_CALIBRATION_ACCESSION # internal calibration
                            calibration_status = :internal
                        elseif acc == INSTRUMENT_SPECIFIC_CALIBRATION_ACCESSION && calibration_status == :uncalibrated # instrument specific calibration
                            calibration_status = :internal # Assume as a form of internal calibration
                        elseif acc == LASER_WAVELENGTH_ACCESSION # laser wavelength
                            laser_settings["wavelength_nm"] = tryparse(Float64, val)
                        elseif acc == LASER_FLUENCE_ACCESSION # laser fluence
                            laser_settings["fluence"] = tryparse(Float64, val)
                        elseif acc == LASER_REPETITION_RATE_ACCESSION # laser repetition rate
                            laser_settings["repetition_rate_hz"] = tryparse(Float64, val)
                        end
                    end
                end
            end
        end
    catch e
        @warn "Could not fully parse instrument metadata from mzML header. Using defaults. Error: $e"
    finally
        seekstart(stream) # Ensure stream is reset after header parsing for subsequent operations
    end
    
    println("DEBUG: Finished mzML instrument metadata parsing.")
    
    # Always return a valid object
    return InstrumentMetadata(
        resolution,
        instrument_meta.acquisition_mode, # To be determined by load_mzml_lazy
        instrument_meta.mz_axis_type,
        calibration_status,
        instrument_model,
        mass_accuracy_ppm,
        isempty(laser_settings) ? nothing : laser_settings,
        polarity,
        isempty(vendor_preprocessing_steps) ? nothing : vendor_preprocessing_steps # NEW field
    )
end

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
    start_pos_bda_block = position(stream) # Mark start of potential BDA block
    
    # Find the <binaryDataArray> tag and extract encodedLength
    # We read the entire tag line into a string, then parse.
    bda_tag_line = find_tag(stream, r"<binaryDataArray\s+encodedLength=\"(\d+)\"") # This is still using `find_tag` for now
    encoded_length = parse(Int32, bda_tag_line.captures[1])
    
    # Mark the position right after the <binaryDataArray> opening tag (and its attributes)
    # This is needed to seek back to find the <binary> tag.
    pos_after_bda_open = position(stream)

    # Now, read the entire block content between the opening <binaryDataArray> and closing </binaryDataArray>
    # This captures all cvParams and the <binary> tag.
    # Note: readuntil consumes the delimiter.
    # To avoid repeated allocations by readline, we read the whole inner block once.
    bda_inner_content = readuntil(stream, "</binaryDataArray>") 
    
    # Initialize parameters
    data_format::DataType = Float64
    compression_flag::Bool = false
    axis::Symbol = :mz

    # Parse bda_inner_content for cvParams
    # Iterate over matches of cvParam tags within the content string
    cvparam_matches = eachmatch(r"<cvParam\s+accession=\"([^\"]+)\"[^>]*>", bda_inner_content)
    for cv_match in cvparam_matches
        acc_str = cv_match.captures[1]
        
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

    # Find binary_offset by seeking back and using readuntil for "<binary>"
    seek(stream, pos_after_bda_open) # Seek to just after the opening <binaryDataArray> tag
    readuntil(stream, "<binary>") # Read until "<binary>"
    binary_offset = position(stream) # Position is now after "<binary>"

    # The stream is already positioned at the end of the </binaryDataArray> block because of the first readuntil.
    # So no need to readuntil again.

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
    
    local id::String = ""
    local mode::SpectrumMode = UNKNOWN
    local mz_asset::Union{SpectrumAsset, Nothing} = nothing
    local int_asset::Union{SpectrumAsset, Nothing} = nothing
    
    # Read line by line to extract information
    spectrum_start_line = readline(stream) # Read the <spectrum> tag
    
    id_match = match(SPECTRUM_TAG_ID_REGEX, spectrum_start_line)
    id = id_match === nothing ? "" : id_match.captures[1]

    while !eof(stream)
        line = readline(stream)

        if occursin("</spectrum>", line)
            break
        end

        # Check for mode
        if mode == UNKNOWN # Only set if not already determined
            accession_match = match(ACCESSION_REGEX, line)
            if accession_match !== nothing
                acc = accession_match.captures[1]
                if acc == CENTROID_SPECTRUM_ACCESSION
                    mode = CENTROID
                elseif acc == PROFILE_SPECTRUM_ACCESSION
                    mode = PROFILE
                end
            end
        end

        # Find binaryDataArrayList and parse assets
        if occursin("<binaryDataArrayList", line)
            # We are now positioned at the line after <binaryDataArrayList>
            # The next two binaryDataArray tags should be immediately following.
            asset1 = get_spectrum_asset_metadata(stream)
            asset2 = get_spectrum_asset_metadata(stream)

            # Determine which asset is m/z and which is intensity
            mz_asset, int_asset = if asset1.axis_type == :mz
                (asset1, asset2)
            else
                (asset2, asset1)
            end
            
            # After parsing binaryDataArrays, we can expect to quickly hit </spectrumList> or other closing tags.
            # We can break early if we have all necessary info, or continue to ensure stream is advanced.
            break # All essential info is found
        end
    end

    if mz_asset === nothing || int_asset === nothing
        # Handle cases where assets weren't found (e.g., empty spectrum or parsing error)
        # Provide default/placeholder assets to avoid errors downstream
        @warn "mzML: Binary assets not found for spectrum $id at offset $offset. Using placeholder assets."
        # These default formats should ideally come from global context if available, otherwise assume typical
        default_mz_format = Float64 
        default_intensity_format = Float64
        mz_asset = SpectrumAsset(default_mz_format, false, Int64(0), 0, :mz)
        int_asset = SpectrumAsset(default_intensity_format, false, Int64(0), 0, :intensity)
    end

    return SpectrumMetadata(Int32(0), Int32(0), id, :sample, mode, mz_asset, int_asset)
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
    
    # Actively search for offset tags, ignoring other lines until the end of the index is found.
    while !eof(stream)
        line = readline(stream)
        
        # First, check for the end condition.
        if occursin("</index>", line) || occursin("</indexedmzML>", line)
            break
        end

        # If it's not the end, see if it's an offset tag.
        m = match(OFFSET_TAG_REGEX, line)
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
    
    index_offset_match = match(INDEX_LIST_OFFSET_REGEX, footer)
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
function load_mzml_lazy(file_path::String; cache_size::Int=100, use_mmap::Bool=false) # use_mmap is ignored for mzML files
    println("DEBUG: Opening file stream for $file_path")
    ts_stream = ThreadSafeFileHandle(file_path, "r")

    try
        # --- NEW: Parse instrument metadata from header ---
        println("DEBUG: Parsing instrument metadata from header...")
        instrument_meta = parse_instrument_metadata_mzml(ts_stream.handle)

        println("--- Extracted Instrument Metadata ---")
        println("Resolution: ", instrument_meta.resolution)
        println("Acquisition Mode (pre-check): ", instrument_meta.acquisition_mode)
        println("Calibration Status: ", instrument_meta.calibration_status)
        println("Instrument Model: ", instrument_meta.instrument_model)
        println("Mass Accuracy (ppm): ", instrument_meta.mass_accuracy_ppm)
        println("Laser Settings: ", instrument_meta.laser_settings)
        println("Polarity: ", instrument_meta.polarity)
        println("------------------------------------")

        seekstart(ts_stream.handle) # Reset stream after header parsing

        println("DEBUG: Finding index offset...")
        index_offset = find_index_offset(ts_stream.handle)
        println("DEBUG: Seeking to index list at offset $index_offset.")
        seek(ts_stream.handle, index_offset)
        
        println("DEBUG: Searching for '<index name=\"spectrum\">'.")
        if find_tag(ts_stream.handle, INDEX_SPECTRUM_TAG_REGEX) === nothing
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
        println("DEBUG: Metadata parsing complete for all $num_spectra spectra.")

        # Assuming uniform data formats, take from the first spectrum
        first_meta = spectra_metadata[1]
        mz_format = first_meta.mz_asset.format
        intensity_format = first_meta.int_asset.format
        println("DEBUG: Inferred global m/z format: $mz_format")
        println("DEBUG: Inferred global intensity format: $intensity_format")
        
        # --- NEW: Determine overall acquisition mode ---
        modes = [meta.mode for meta in spectra_metadata]
        num_centroid = count(m -> m == CENTROID, modes)
        num_profile = count(m -> m == PROFILE, modes)
        
        acq_mode_symbol = if num_centroid > 0 && num_profile == 0
            :centroid
        elseif num_profile > 0 && num_centroid == 0
            :profile
        elseif num_centroid > 0 && num_profile > 0
            :mixed
        else
            :unknown
        end
        println("DEBUG: Inferred overall acquisition mode: $acq_mode_symbol (Centroid: $num_centroid, Profile: $num_profile)")

        final_instrument_meta = InstrumentMetadata(
            instrument_meta.resolution,
            acq_mode_symbol, # Update with parsed mode
            instrument_meta.mz_axis_type,
            instrument_meta.calibration_status,
            instrument_meta.instrument_model,
            instrument_meta.mass_accuracy_ppm,
            instrument_meta.laser_settings,
            instrument_meta.polarity,
            instrument_meta.vendor_preprocessing_steps # Add this new field
        )

        source = MzMLSource(ts_stream, mz_format, intensity_format)
        println("DEBUG: Creating MSIData object.")
        return MSIData(source, spectra_metadata, final_instrument_meta, (0, 0), nothing, cache_size)

    catch e
        close(ts_stream) # Ensure stream is closed on error
        rethrow(e)
    end
end
