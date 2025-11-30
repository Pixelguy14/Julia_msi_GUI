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
        # Reset stream and read a sufficiently large header block to find metadata.
        seekstart(stream)
        header_block = ""
        header_read_limit = 50000 # Read up to 50KB of header
        bytes_read = 0
        
        while !eof(stream) && bytes_read < header_read_limit
            line_pos = position(stream)
            line = readline(stream)
            bytes_read += (position(stream) - line_pos)
            
            # Stop at the start of the main data section
            if occursin("<run>", line) || occursin("<spectrumList>", line)
                break
            end
            header_block *= line * "\n"
        end

        header_io = IOBuffer(header_block)
        while !eof(header_io)
            line = readline(header_io)
            
            if occursin("<cvParam", line)
                accession_match = get_attribute(line, "accession")
                value_match = get_attribute(line, "value")
                name_match = get_attribute(line, "name") # For laser attributes and vendor_preprocessing names
                
                if accession_match !== nothing
                    acc = accession_match.captures[1]
                    val = (value_match !== nothing) ? value_match.captures[1] : ""
                    name = (name_match !== nothing) ? name_match.captures[1] : "" # NEW: Get name for logging

                    if acc == "MS:1000031" # instrument model
                        instrument_model = val
                        #println("DEBUG: Instrument model: $instrument_model")
                    elseif acc == "MS:1001496" # mass resolving power (more specific)
                        resolution = tryparse(Float64, val)
                        #println("DEBUG: Resolution (specific): $resolution")
                    elseif acc == "MS:1000011" && resolution === nothing # resolution (less specific)
                        resolution = tryparse(Float64, val)
                        #println("DEBUG: Resolution (less specific): $resolution")
                    elseif acc == "MS:1000016" # mass accuracy (ppm)
                        mass_accuracy_ppm = tryparse(Float64, val)
                        #println("DEBUG: Mass accuracy (ppm): $mass_accuracy_ppm")
                    elseif acc == "MS:1000130" # positive scan
                        polarity = :positive
                        #println("DEBUG: Polarity: positive")
                    elseif acc == "MS:1000129" # negative scan
                        polarity = :negative
                        #println("DEBUG: Polarity: negative")
                    elseif acc == "MS:1000592" # external calibration
                        calibration_status = :external
                        #println("DEBUG: Calibration status: external")
                    elseif acc == "MS:1000593" # internal calibration
                        calibration_status = :internal
                        #println("DEBUG: Calibration status: internal")
                    elseif acc == "MS:1000747" && calibration_status == :uncalibrated # instrument specific calibration
                        calibration_status = :internal # Assume as a form of internal calibration
                        #println("DEBUG: Calibration status: instrument specific (internal)")
                    elseif acc == "MS:1000867" # laser wavelength
                        laser_settings["wavelength_nm"] = tryparse(Float64, val)
                        #println("DEBUG: Laser wavelength: $(laser_settings["wavelength_nm"]) nm")
                    elseif acc == "MS:1000868" # laser fluence
                        laser_settings["fluence"] = tryparse(Float64, val)
                        #println("DEBUG: Laser fluence: $(laser_settings["fluence"])")
                    elseif acc == "MS:1000869" # laser repetition rate
                        laser_settings["repetition_rate_hz"] = tryparse(Float64, val)
                        #println("DEBUG: Laser repetition rate: $(laser_settings["repetition_rate_hz"]) Hz")
                    # NEW: Vendor Preprocessing terms
                    elseif acc == "MS:1000579" # baseline correction
                        push!(vendor_preprocessing_steps, "Baseline Correction")
                        #println("DEBUG: Vendor preprocessing step: Baseline Correction")
                    elseif acc == "MS:1000580" # smoothing
                        push!(vendor_preprocessing_steps, "Smoothing")
                        #println("DEBUG: Vendor preprocessing step: Smoothing")
                    elseif acc == "MS:1000578" # data transformation (e.g., centroiding)
                        push!(vendor_preprocessing_steps, "Data Transformation: $(name)")
                        #println("DEBUG: Vendor preprocessing step: Data Transformation: $(name)")
                    elseif acc == "MS:1000800" # deisotoping
                        push!(vendor_preprocessing_steps, "Deisotoping")
                        #println("DEBUG: Vendor preprocessing step: Deisotoping")
                    end
                end
            end
        end
    catch e
        @warn "Could not fully parse instrument metadata from mzML header. Using defaults. Error: $e"
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
    start_pos = position(stream)
    
    #println("DEBUG: Entering get_spectrum_asset_metadata to parse binaryDataArray...")

    bda_tag = find_tag(stream, r"<binaryDataArray\s+encodedLength=\"(\d+)\"")

    if bda_tag === nothing
        throw(FileFormatError("Cannot find binaryDataArray"))
    end
    encoded_length = parse(Int32, bda_tag.captures[1])
    #println("DEBUG:   Encoded length: $encoded_length")

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
                #println("DEBUG:   Axis type identified as: m/z")
            elseif acc_str == INTENSITY_AXIS_ACCESSION
                axis = :intensity
                #println("DEBUG:   Axis type identified as: intensity")
            elseif haskey(DATA_FORMAT_ACCESSIONS, acc_str)
                data_format = DATA_FORMAT_ACCESSIONS[acc_str]
                #println("DEBUG:   Data format identified as: $data_format")
            elseif acc_str == COMPRESSION_ACCESSION
                compression_flag = true
                #println("DEBUG:   Compression: true")
            elseif acc_str == NO_COMPRESSION_ACCESSION
                compression_flag = false
                #println("DEBUG:   Compression: false")
            end
        end
    end

    seek(stream, start_pos)
    readuntil(stream, "<binary>")
    binary_offset = position(stream)
    #println("DEBUG:   Binary data offset: $binary_offset")

    # Move stream to the end of the binary data array for the next iteration
    readuntil(stream, "</binaryDataArray>")
    #println("DEBUG: Exiting get_spectrum_asset_metadata.")

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
    #println("DEBUG:   Parsing spectrum ID: $id")

    # Determine mode from the XML block
    mode = UNKNOWN
    if occursin("MS:1000127", spectrum_xml)
        mode = CENTROID
        #println("DEBUG:   Spectrum mode: CENTROID")
    elseif occursin("MS:1000128", spectrum_xml)
        mode = PROFILE
        #println("DEBUG:   Spectrum mode: PROFILE")
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

    #println("DEBUG:     m/z Asset - Format: $(mz_asset.format), Compressed: $(mz_asset.is_compressed), Offset: $(mz_asset.offset), Encoded Length: $(mz_asset.encoded_length)")
    #println("DEBUG:     Intensity Asset - Format: $(int_asset.format), Compressed: $(int_asset.is_compressed), Offset: $(int_asset.offset), Encoded Length: $(int_asset.encoded_length)")
    #println("DEBUG:   Finished parsing spectrum ID: $id metadata.")

    # Create the new unified metadata object
    # For mzML, x and y coordinates are not applicable, so we use 0.
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
