"""
This file contains the workflow for converting .mzML files (with one spectrum per pixel)
into a proper .imzML/.ibd file pair, using a separate synchronization file.
It replicates the functionality of the original R scripts that use MALDIquant.
"""

using DataFrames, Printf, CSV

# This file assumes that the main application file (e.g., app.jl) has already included
# the necessary source files: MSIData.jl, mzML.jl, imzML.jl


# A struct to hold the processed pixel data before exporting
struct ProcessedPixel
    coords::Tuple{Int, Int}
    mz::Vector{Float64} # Assuming m/z is consistent, can be optimized later
    intensity::Vector{Float32}
end

# This struct will hold the metadata needed for the XML file
struct BinaryMetadata
    mz_offset::UInt64
    mz_length::UInt64
    int_offset::UInt64
    int_length::UInt64
end

"""
    GetMzmlScanTime(fileName::String)

Parses a .mzML file to extract the scan start time for each spectrum.

# Arguments
* `fileName`: Path to the .mzML file.

# Returns
- A `Matrix{Int64}` where each row is `[spectrum_index, time_in_milliseconds]`.
"""
function GetMzmlScanTime(fileName::String)
    times = Tuple{Int64, Int64}[]
    
    # Pre-allocate based on an estimate to reduce re-allocations
    try
        file_size = filesize(fileName)
        estimated_spectra = max(1000, file_size ÷ 10000) # Heuristic
        sizehint!(times, estimated_spectra)
    catch
        # Ignore if filesize fails, sizehint! is just an optimization
    end

    open(fileName, "r") do stream
        state = :outside_spectrum
        current_index = 0
        current_time = nothing
        
        for line in eachline(stream)
            if state == :outside_spectrum
                if occursin("<spectrum ", line)
                    state = :in_spectrum
                    current_index = 0
                    current_time = nothing
                    
                    # Process index from the start tag itself
                    idx_match = match(r"index=\"(\d+)\"", line)
                    if idx_match !== nothing
                        current_index = parse(Int, idx_match.captures[1]) + 1
                    end
                end
            elseif state == :in_spectrum
                # Look for time cvParam
                if occursin("<cvParam", line)
                    if occursin("MS:1000016", line) || occursin("MS:1000015", line)
                        time_match = match(r"value=\"([\d\.]+)\"", line)
                        if time_match !== nothing
                            time_val = parse(Float64, time_match.captures[1])
                            # MS:1000016 is minutes, MS:1000015 is seconds
                            unit_scale = occursin("MS:1000016", line) ? 60000 : 1000
                            current_time = round(Int64, time_val * unit_scale)
                        end
                    end
                end

                # Check for end of spectrum
                if occursin("</spectrum>", line)
                    state = :outside_spectrum
                    if current_index > 0 && current_time !== nothing
                        push!(times, (current_index, current_time))
                    end
                end
            end
        end
    end
    
    # Sort by spectrum index to ensure correct order before normalization
    sort!(times, by = x -> x[1])

    # Normalize times relative to the first scan
    if !isempty(times)
        first_time = times[1][2]
        for i in eachindex(times)
            times[i] = (times[i][1], times[i][2] - first_time)
        end
    end
    
    # Convert vector of tuples to a matrix
    if isempty(times)
        return Matrix{Int64}(undef, 0, 2)
    end
    return permutedims(hcat(collect.(times)...))
end

"""
    MatchAcquireTime(sync_file_path::String, scans::Matrix{Int64})

Correlates pixel acquisition times from a synchronization file with scan
acquisition times from an mzML file.

This is a Julia implementation of the `MatchAcquireTime` function from the R scripts.

# Arguments
* `sync_file_path`: Path to the synchronization file (.txt).
* `scans`: A matrix of scan times, as returned by `GetMzmlScanTime`.

# Returns
- A `Matrix{Int64}` where each row is `[pixel_index, pixel_time_ms, first_scan_index, last_scan_index]`.
"""
function MatchAcquireTime(sync_file_path::String, scans::Matrix{Int64}; img_width::Int=0, img_height::Int=0)
    if !isfile(sync_file_path)
        error("Synchronization file not found: $sync_file_path")
    end

    pixel_df = CSV.read(sync_file_path, DataFrame, header=false, skipto=3)
    pixel_matrix = Matrix(pixel_df)
    
    num_pixels_original = size(pixel_matrix, 1)
    if num_pixels_original == 0
        return zeros(Int64, 0, 5)
    end

    # Determine image dimensions and generate coordinates
    current_img_width = img_width
    current_img_height = img_height

    if current_img_width == 0 || current_img_height == 0
        if size(pixel_matrix, 2) >= 3
            @info "Sync file contains X, Y coordinates"
            coordinates = convert(Matrix{Int}, pixel_matrix[:, 1:2])
            pixel_times = convert(Vector{Int64}, pixel_matrix[:, 3])
            current_img_width = maximum(coordinates[:, 1])
            current_img_height = maximum(coordinates[:, 2])
            
            # Apply R's pixel truncation logic
            num_pixels = (floor(Int, (num_pixels_original - 1) / current_img_width)) * current_img_width
            if num_pixels <= 0
                error("Calculated num_pixels to process is zero or negative: $num_pixels")
            end
            
            # Truncate arrays
            final_coordinates = coordinates[1:num_pixels, :]
            final_pixel_times = pixel_times[1:num_pixels]
            
        else
            @info "Sync file contains Index and Time; generating coordinates"
            pixel_times = convert(Vector{Int64}, pixel_matrix[:, 2])
            
            # R's width detection logic
            diffs = pixel_matrix[2:end, 1] .- pixel_matrix[1:end-1, 1]
            width_indices = findall(x -> x != 1, diffs)
            if !isempty(width_indices)
                current_img_width = width_indices[1]
            else
                current_img_width = num_pixels_original
            end
            current_img_height = num_pixels_original ÷ current_img_width
            
            # R's pixel truncation logic
            num_pixels = (floor(Int, (num_pixels_original - 1) / current_img_width)) * current_img_width
            if num_pixels <= 0
                error("Calculated num_pixels to process is zero or negative: $num_pixels")
            end
            
            # Generate coordinates for truncated pixels
            final_coordinates = zeros(Int, num_pixels, 2)
            final_pixel_times = zeros(Int64, num_pixels)
            
            for i in 1:num_pixels
                idx = pixel_matrix[i, 1]
                final_coordinates[i, 1] = ((idx - 1) % current_img_width) + 1
                final_coordinates[i, 2] = fld(idx - 1, current_img_width) + 1
                final_pixel_times[i] = pixel_times[i]
            end
        end
    else
        @info "Using provided dimensions: $(current_img_width)x$(current_img_height)"
        pixel_times = convert(Vector{Int64}, pixel_matrix[:, 2])
        
        # R's pixel truncation logic (FIXED: consistent formula)
        num_pixels = (floor(Int, (num_pixels_original - 1) / current_img_width)) * current_img_width
        if num_pixels <= 0
            error("Calculated num_pixels to process is zero or negative: $num_pixels")
        end
        
        # Generate coordinates
        final_coordinates = zeros(Int, num_pixels, 2)
        final_pixel_times = zeros(Int64, num_pixels)
        
        for i in 1:num_pixels
            idx = pixel_matrix[i, 1]
            final_coordinates[i, 1] = ((idx - 1) % current_img_width) + 1
            final_coordinates[i, 2] = fld(idx - 1, current_img_width) + 1
            final_pixel_times[i] = pixel_times[i]
        end
    end

    # Normalize pixel times
    if !isempty(final_pixel_times)
        min_pixel_time = minimum(final_pixel_times)
        final_pixel_times .-= min_pixel_time
    end

    num_scans = size(scans, 1)
    if num_scans == 0
        return hcat(final_coordinates, final_pixel_times, zeros(Int64, num_pixels, 2))
    end
    
    # FIXED R's time matching algorithm with bounds checking
    first_idx = 1
    last_idx = 1
    index_matrix = zeros(Int, num_pixels, 2)
    
    for i_pixel in 1:num_pixels
        pixel_time = final_pixel_times[i_pixel]
        
        # Find the first scan that reaches or exceeds pixel time
        while last_idx <= num_scans && scans[last_idx, 2] < pixel_time
            last_idx += 1
        end
        
        # Ensure valid indices
        if last_idx > num_scans
            # No more scans available for remaining pixels
            index_matrix[i_pixel:end, 1] .= num_scans + 1  # Invalid index
            index_matrix[i_pixel:end, 2] .= num_scans      # Invalid index
            break
        end
        
        # Assign scan indices with bounds checking
        start_scan = max(1, first_idx)
        end_scan = max(1, last_idx - 1)
        
        # Ensure start_scan <= end_scan
        if start_scan > end_scan
            start_scan = end_scan
        end
        
        index_matrix[i_pixel, 1] = start_scan
        index_matrix[i_pixel, 2] = end_scan
        
        # Update for next iteration (R's algorithm)
        first_idx = last_idx - 1
        last_idx = first_idx
        
        # Ensure first_idx doesn't go below 1
        if first_idx < 1
            first_idx = 1
            last_idx = 1
        end
    end
    
    return hcat(final_coordinates, final_pixel_times, index_matrix)
end


function RenderPixel(pixel_info, scans, msi_data::MSIData, scan_time_deltas, pixel_time_deltas)
    pixel_time = pixel_info[3]
    first_scan = pixel_info[4]
    last_scan = pixel_info[5]
    num_actions = last_scan - first_scan

    # Get reference m/z array from the first scan involved
    mz_array, _ = GetSpectrum(msi_data, first_scan)
    new_intensity = zeros(Float32, length(mz_array))

    # SAFETY: Ensure we have valid scan indices
    if first_scan < 1 || last_scan > size(scans, 1) || first_scan > last_scan
        return (mz_array, new_intensity)  # Return zero intensity for invalid ranges
    end

    if num_actions == 0 # Single scan contributes to the pixel
        _, intensity = GetSpectrum(msi_data, first_scan)
        # SAFETY: Ensure positive scaling
        scale = max(pixel_time_deltas[first_scan] / scan_time_deltas[first_scan], 0.0f0)
        new_intensity .= intensity .* scale

    elseif num_actions == 1 # Two partial scans contribute
        # First partial scan
        _, intensity1 = GetSpectrum(msi_data, first_scan)
        scale1 = max((scans[last_scan, 2] - pixel_time) / scan_time_deltas[first_scan], 0.0f0)
        new_intensity .+= intensity1 .* scale1

        # Second partial scan
        _, intensity2 = GetSpectrum(msi_data, last_scan)
        next_pixel_time = pixel_time + pixel_time_deltas[first_scan]
        scale2 = max((next_pixel_time - scans[last_scan, 2]) / scan_time_deltas[last_scan], 0.0f0)
        new_intensity .+= intensity2 .* scale2

    elseif num_actions > 1 # Multiple scans contribute
        # First partial scan
        _, intensity1 = GetSpectrum(msi_data, first_scan)
        scale1 = max((scans[first_scan + 1, 2] - pixel_time) / scan_time_deltas[first_scan], 0.0f0)
        new_intensity .+= intensity1 .* scale1

        # Full scans in the middle
        for i in (first_scan + 1):(last_scan - 1)
            _, intensity_middle = GetSpectrum(msi_data, i)
            new_intensity .+= intensity_middle
        end

        # Last partial scan
        _, intensity2 = GetSpectrum(msi_data, last_scan)
        next_pixel_time = pixel_time + pixel_time_deltas[first_scan]
        scale2 = max((next_pixel_time - scans[last_scan, 2]) / scan_time_deltas[last_scan], 0.0f0)
        new_intensity .+= intensity2 .* scale2
    end

    # FINAL SAFETY: Clamp any negative values to zero
    new_intensity = max.(new_intensity, 0.0f0)
    
    return (mz_array, new_intensity)
end

function ConvertMzmlToImzml(source_file::String, timing_matrix::Matrix{Int64}, scans::Matrix{Int64})
    if size(timing_matrix, 1) == 0
        return ProcessedPixel[], (0, 0)
    end
    
    # Get image dimensions from timing_matrix (assuming columns 1 and 2 are X and Y)
    width = maximum(timing_matrix[:, 1])
    height = maximum(timing_matrix[:, 2])

    # Open the mzML file using the new memory-efficient API
    msi_data = OpenMSIData(source_file)

    # Pre-calculate time deltas robustly
    scan_time_deltas = zeros(Int64, size(scans, 1))
    if size(scans, 1) > 1
        for i in 1:(size(scans, 1) - 1)
            delta = scans[i+1, 2] - scans[i, 2]
            scan_time_deltas[i] = max(1, delta) # Ensure delta is at least 1
        end
        scan_time_deltas[end] = max(1, scan_time_deltas[end-1]) # Use previous delta for last scan, ensure at least 1
    end

    pixel_time_deltas = zeros(Int64, size(timing_matrix, 1))
    if size(timing_matrix, 1) > 1
        for i in 1:(size(timing_matrix, 1) - 1)
            # Assumes time is in the 3rd column of the timing_matrix
            delta = timing_matrix[i+1, 3] - timing_matrix[i, 3]
            pixel_time_deltas[i] = max(1, delta) # Ensure delta is at least 1
        end
        pixel_time_deltas[end] = max(1, pixel_time_deltas[end-1]) # Use previous delta for last pixel, ensure at least 1
    end

    processed_pixels = ProcessedPixel[]
    sizehint!(processed_pixels, size(timing_matrix, 1))
    empty_pixel_count = 0

    # DEBUG: Scan coverage analysis
    println("DEBUG: Scan coverage analysis:")
    covered_pixels_count = 0
    for i in 1:size(timing_matrix, 1)
        first_scan = timing_matrix[i, 4]
        last_scan = timing_matrix[i, 5]
        if first_scan <= last_scan && first_scan >= 1 && last_scan <= size(scans, 1)
            covered_pixels_count += 1
        end
    end
    println("  Pixels with potential scans: $covered_pixels_count/$(size(timing_matrix, 1))")
    println("  Total scans available: $(size(scans, 1))")

    for i in 1:size(timing_matrix, 1)
        pixel_info = timing_matrix[i, :]
        x = pixel_info[1]
        y = pixel_info[2]
        first_scan = pixel_info[4]
        last_scan = pixel_info[5]

        if first_scan > last_scan || first_scan < 1 || last_scan > size(scans, 1)
            empty_pixel_count += 1
            push!(processed_pixels, ProcessedPixel((x, y), Float64[], Float32[]))
            continue
        end

        mz, intensity = RenderPixel(pixel_info, scans, msi_data, scan_time_deltas, pixel_time_deltas)
        push!(processed_pixels, ProcessedPixel((x, y), mz, intensity))
    end

    @info "Found and created $empty_pixel_count empty pixels out of $(size(timing_matrix, 1)) total."
    return processed_pixels, (width, height)
end

function ExportImzml(target_file::String, pixels::Vector{ProcessedPixel}, dims::Tuple{Int, Int})
    ibd_file = replace(target_file, r"\.imzML$"i => ".ibd")
    
    # Ensure we have valid pixels
    if isempty(pixels)
        @error "No pixels to export"
        return false
    end

    binary_meta = BinaryMetadata[]
    sizehint!(binary_meta, length(pixels))

    try
        # Step 1: Write .ibd file with proper binary format
        open(ibd_file, "w") do ibd_stream
            # Write UUID as first 16 bytes (placeholder)
            write(ibd_stream, zeros(UInt8, 16))

            # Write in little-endian format (standard for imzML)
            for pixel in pixels
                # Write m/z array as 64-bit little-endian floats
                mz_offset = position(ibd_stream)
                for val in pixel.mz
                    write(ibd_stream, htol(Float64(val)))  # Host to little-endian
                end
                mz_length = position(ibd_stream) - mz_offset

                # Write intensity array as 32-bit little-endian floats  
                int_offset = position(ibd_stream)
                for val in pixel.intensity
                    write(ibd_stream, htol(Float32(val)))  # Host to little-endian
                end
                int_length = position(ibd_stream) - int_offset
                
                push!(binary_meta, BinaryMetadata(mz_offset, mz_length, int_offset, int_length))
            end
        end

        # Step 2: Write proper .imzML XML file
        open(target_file, "w") do imzml_stream
            # Use indexedmzML wrapper
            write(imzml_stream, """<?xml version="1.0" encoding="ISO-8859-1"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd">
  <mzML version="1.1" id="$(splitext(basename(target_file))[1])">
""")

            # CV List
            write(imzml_stream, """  <cvList count="3">
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="1.3.1" URI="http://psidev.info/ms/mzML/psi-ms.obo"/>
    <cv id="UO" fullName="Unit Ontology" version="1.15" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>
    <cv id="IMS" fullName="Imaging MS Ontology" version="0.9.1" URI="http://www.maldi-msi.org/download/imzml/imagingMS.obo"/>
  </cvList>
""")

            # File Description
            write(imzml_stream, """  <fileDescription>
    <fileContent>
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="mass spectrum"/>
      <cvParam cvRef="IMS" accession="IMS:1000031" name="processed"/>
    </fileContent>
  </fileDescription>
""")

            # Referenceable Param Groups
            write(imzml_stream, """  <referenceableParamGroupList count="2">
    <referenceableParamGroup id="mzArray">
      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensityArray">
      <cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
      <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/>
    </referenceableParamGroup>
  </referenceableParamGroupList>
""")

            # Sample List
            write(imzml_stream, """  <sampleList count="1">
    <sample id="sample1" name="ImagingSample">
      <cvParam cvRef="MS" accession="MS:1000001" name="sample number" value="1"/>
    </sample>
  </sampleList>
""")

            # Software List - OPEN SOURCE
            write(imzml_stream, """  <softwareList count="1">
    <software id="MSIConverter" version="1.0">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="mzML to imzML converter"/>
    </software>
  </softwareList>
""")

            # Scan Settings IMAGING
            write(imzml_stream, """  <scanSettingsList count="1">
    <scanSettings id="scanSettings1">
      <cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixel x" value="$(dims[1])"/>
      <cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixel y" value="$(dims[2])"/>
      <cvParam cvRef="IMS" accession="IMS:1000046" name="pixel size x" value="1" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
      <cvParam cvRef="IMS" accession="IMS:1000047" name="pixel size y" value="1" unitCvRef="UO" unitAccession="UO:0000017" unitName="micrometer"/>
    </scanSettings>
  </scanSettingsList>
""")

            # Instrument Configuration MALDI/TOF
            write(imzml_stream, """  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="instrument1">
      <componentList count="3">
        <source order="1">
          <cvParam cvRef="MS" accession="MS:1000075" name="MALDI source"/>
        </source>
        <analyzer order="2">
          <cvParam cvRef="MS" accession="MS:1000084" name="time-of-flight"/>
        </analyzer>
        <detector order="3">
            <cvParam cvRef="MS" accession="MS:1000253" name="electron multiplier"/>
        </detector>
      </componentList>
      <softwareRef ref="MSIConverter"/>
    </instrumentConfiguration>
  </instrumentConfigurationList>
""")

            # Data Processing
            write(imzml_stream, """  <dataProcessingList count="1">
    <dataProcessing id="conversionProcessing">
      <processingMethod order="1" softwareRef="MSIConverter">
        <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to imzML"/>
      </processingMethod>
    </dataProcessing>
  </dataProcessingList>
""")

            # Run and Spectrum List
            spectrum_offsets = UInt64[]
            write(imzml_stream, """  <run defaultInstrumentConfigurationRef="instrument1" id="run1" sampleRef="sample1">
    <spectrumList count="$(length(pixels))" defaultDataProcessingRef="conversionProcessing">
""")

            # Write each spectrum
            for (i, pixel) in enumerate(pixels)
                meta = binary_meta[i]
                x, y = pixel.coords
                
                spectrum_start = position(imzml_stream)
                push!(spectrum_offsets, spectrum_start)

                write(imzml_stream, """      <spectrum id="Scan=$(i)" defaultArrayLength="$(length(pixel.mz))" index="$(i-1)">
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum"/>
        <scanList count="1">
          <scan instrumentConfigurationRef="instrument1">
            <cvParam cvRef="IMS" accession="IMS:1000050" name="position x" value="$(x)"/>
            <cvParam cvRef="IMS" accession="IMS:1000051" name="position y" value="$(y)"/>
          </scan>
        </scanList>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mzArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="$(meta.mz_offset)"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="$(length(pixel.mz))"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="$(meta.mz_length)"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="intensityArray"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external offset" value="$(meta.int_offset)"/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external array length" value="$(length(pixel.intensity))"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external encoded length" value="$(meta.int_length)"/>
            <binary/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
""")
            end

            write(imzml_stream, """    </spectrumList>
  </run>
</mzML>
""")
            # Index List
            index_list_start = position(imzml_stream)
            write(imzml_stream, """  <indexList count="1">
    <index name="spectrum">
""")
            for (i, offset) in enumerate(spectrum_offsets)
                write(imzml_stream, "      <offset idRef=\"Scan=$(i)\">$offset</offset>\n")
            end
            write(imzml_stream, """    </index>
  </indexList>
  <indexListOffset>$index_list_start</indexListOffset>
</indexedmzML>
""")
        end

        println("Successfully created: $target_file")
        println("Successfully created: $ibd_file")
        return true
        
    catch e
        @error "Failed to export imzML files" exception=(e, catch_backtrace())
        # Clean up partial files
        isfile(ibd_file) && rm(ibd_file, force=true)
        isfile(target_file) && rm(target_file, force=true)
        return false
    end
end

"""
    ImportMzmlFile(source_file::String, sync_file::String, target_file::String)

Main workflow function to convert a .mzML file to an .imzML file.

# Arguments
* `source_file`: Path to the input .mzML file.
* `sync_file`: Path to the synchronization text file.
* `target_file`: Path for the output .imzML file (the .ibd will be named accordingly).
* `img_width`: width dimention for the creation of the x axis
* `img_height`: height dimention for the creation of the y axis
"""
function ImportMzmlFile(source_file::String, sync_file::String, target_file::String; img_width::Int=0, img_height::Int=0)
    println("Step 1: Getting scan times from .mzML file...")
    scans = GetMzmlScanTime(source_file)

    println("Step 2: Matching acquisition times...")
    timing_matrix = MatchAcquireTime(sync_file, scans; img_width=img_width, img_height=img_height)

    println("Step 3: Converting spectra...")
    processed_pixels, (width, height) = ConvertMzmlToImzml(source_file, timing_matrix, scans)

    # Flip image vertically to match R script output
    for i in eachindex(processed_pixels)
        x, y = processed_pixels[i].coords
        processed_pixels[i] = ProcessedPixel((x, height - y + 1), processed_pixels[i].mz, processed_pixels[i].intensity)
    end

    println("Step 4: Exporting to .imzML/.ibd format...")
    success = ExportImzml(target_file, processed_pixels, (width, height))

    if success
        println("Conversion successful: $target_file")
    else
        println("Conversion failed.")
    end
    return success
end
