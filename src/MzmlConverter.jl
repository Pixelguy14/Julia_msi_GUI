# src/MzmlConverter.jl

"""
This file contains the workflow for converting .mzML files (with one spectrum per pixel)
into a proper .imzML/.ibd file pair, using a separate synchronization file.
It replicates the functionality of the project's original R scripts that use MALDIquant.
"""

using CSV, DataFrames

# Note: This file assumes that `LoadSpectra` from DataAccess.jl and the necessary
# structs/functions from ParserHelpers.jl, mzML.jl, and imzML.jl are available
# in the execution context (e.g., included in the main app.jl).

# A struct to hold the processed pixel data before exporting
struct ProcessedPixel
    coords::Tuple{Int, Int}
    mz::Vector{Float64} # Assuming m/z is consistent, can be optimized later
    intensity::Vector{Float32}
end

"""
    GetMzmlScanTime(fileName::String)

Parses a .mzML file to extract the scan start time for each spectrum.
This is a Julia implementation of the `GetMzmlScanTime` function from the R scripts.

# Arguments
* `fileName`: Path to the .mzML file.

# Returns
- A `Matrix{Int64}` where each row is `[spectrum_index, time_in_milliseconds]`.
"""
function GetMzmlScanTime(fileName::String)
    spec_count = 0
    times = Tuple{Int64, Int64}[]

    open(fileName, "r") do stream
        file_content = read(stream, String)
        
        # First, find the total number of spectra to pre-allocate
        spectrum_list_match = match(r"<spectrumList count=\"(\d+)\"", file_content)
        if spectrum_list_match === nothing
            error("Could not find <spectrumList> tag or count attribute.")
        end
        total_spectra = parse(Int, spectrum_list_match.captures[1])
        sizehint!(times, total_spectra)

        # Iterate over each spectrum block
        for spectrum_match in eachmatch(r"<spectrum index=\"(\d+)\"[^>]*>.*?<\/spectrum>", file_content, overlay=true)
            spec_block = spectrum_match.match
            # The R script appears to treat the 0-based mzML index as 1-based, so we add 1.
            spec_index = parse(Int, spectrum_match.captures[1]) + 1

            # Find the scan start time within the spectrum block
            time_match = match(r"<cvParam accession=\"MS:1000016\".*?value=\"([\d\.]+)\"", spec_block)
            
            if time_match !== nothing
                # The value is in minutes. Convert to milliseconds and round.
                time_in_minutes = parse(Float64, time_match.captures[1])
                time_in_ms = round(Int64, time_in_minutes * 60000)
                push!(times, (spec_index, time_in_ms))
            end
        end
        
        # Normalize times so the first scan is at t=0
        if !isempty(times)
            sort!(times, by = x -> x[1]) # Ensure scans are sorted by index
            first_time = times[1][2]
            for i in eachindex(times)
                times[i] = (times[i][1], times[i][2] - first_time)
            end
        end
    end # open

    # Convert vector of tuples to a matrix
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
function MatchAcquireTime(sync_file_path::String, scans::Matrix{Int64})
    if !isfile(sync_file_path)
        error("Synchronization file not found: $sync_file_path")
    end

    # Load xyz stage pixel delays from the sync file, skipping the first 2 header lines.
    pixel_df = CSV.read(sync_file_path, DataFrame, header=false, skipto=3)
    pixel = Matrix(pixel_df)
    
    # Normalize pixel times so the first is at t=0
    pixel[:, 2] .-= pixel[1, 2]

    num_pixels = size(pixel, 1)
    num_scans = size(scans, 1)
    
    index_matrix = zeros(Int64, num_pixels, 2)
    
    scan_row = 1
    for pixel_row in 1:(num_pixels - 1)
        start_scan = scan_row
        
        # Find the scan that ends after the *next* pixel starts
        while scan_row <= num_scans && scans[scan_row, 2] < pixel[pixel_row + 1, 2]
            scan_row += 1
        end
        
        index_matrix[pixel_row, 1] = start_scan
        index_matrix[pixel_row, 2] = scan_row - 1
        
        # The next pixel's scans start from the one that crossed the boundary
        scan_row = max(1, scan_row - 1)
    end
    
    # Assign scans for the last pixel
    index_matrix[num_pixels, 1] = scan_row
    index_matrix[num_pixels, 2] = num_scans
    
    # Handle cases where some pixels have no scans by creating an empty range
    for i in 1:num_pixels
        if index_matrix[i, 2] < index_matrix[i, 1]
            index_matrix[i, 2] = index_matrix[i, 1] - 1
        end
    end

    return hcat(pixel, index_matrix)
end


function RenderPixel(pixel_info, scans, spectra, scan_time_deltas, pixel_time_deltas)
    first_scan = pixel_info[3]
    last_scan = pixel_info[4]
    num_actions = last_scan - first_scan

    # Use the m/z array of the first scan as the reference
    mz_array = spectra[1, first_scan]
    new_intensity = zeros(Float32, length(mz_array))

    if num_actions == 0 # Single scan contributes to the pixel
        scale = pixel_time_deltas[pixel_info[1]] / scan_time_deltas[first_scan]
        new_intensity .+= spectra[2, first_scan] .* scale

    elseif num_actions == 1 # Two partial scans contribute
        # First partial scan
        scale1 = (scans[last_scan, 2] - pixel_info[2]) / scan_time_deltas[first_scan]
        new_intensity .+= spectra[2, first_scan] .* scale1

        # Second partial scan
        next_pixel_time = pixel_info[2] + pixel_time_deltas[pixel_info[1]]
        scale2 = (next_pixel_time - scans[last_scan, 2]) / scan_time_deltas[last_scan]
        new_intensity .+= spectra[2, last_scan] .* scale2

    elseif num_actions > 1 # Multiple scans contribute
        # First partial scan
        scale1 = (scans[first_scan + 1, 2] - pixel_info[2]) / scan_time_deltas[first_scan]
        new_intensity .+= spectra[2, first_scan] .* scale1

        # Full scans in the middle
        for i in (first_scan + 1):(last_scan - 1)
            new_intensity .+= spectra[2, i]
        end

        # Last partial scan
        next_pixel_time = pixel_info[2] + pixel_time_deltas[pixel_info[1]]
        scale2 = (next_pixel_time - scans[last_scan, 2]) / scan_time_deltas[last_scan]
        new_intensity .+= spectra[2, last_scan] .* scale2
    end

    return (mz_array, new_intensity)
end

function ConvertMzmlToImzml(source_file::String, timing_matrix::Matrix{Int64}, scans::Matrix{Int64})
    # Get image width from timing info
    width = findfirst(i -> timing_matrix[i, 1] - timing_matrix[i-1, 1] != 1, 2:size(timing_matrix, 1))
    height = size(timing_matrix, 1) ÷ width

    # Load the full mzML data
    spectra = LoadSpectra(source_file)

    # Pre-calculate time deltas
    scan_time_deltas = diff(scans[:, 2])
    pixel_time_deltas = diff(timing_matrix[:, 2])
    # Append a final delta for the last element
    push!(scan_time_deltas, scan_time_deltas[end])
    push!(pixel_time_deltas, pixel_time_deltas[end])

    processed_pixels = ProcessedPixel[]
    sizehint!(processed_pixels, size(timing_matrix, 1))

    for i in 1:size(timing_matrix, 1)
        pixel_info = timing_matrix[i, :]
        first_scan = pixel_info[3]
        last_scan = pixel_info[4]

        # Skip if there are no scans for this pixel
        (first_scan > last_scan) && continue

        # Calculate coordinates
        x = ((pixel_info[1] - 1) % width) + 1
        y = fld(pixel_info[1] - 1, width) + 1

        mz, intensity = RenderPixel(pixel_info, scans, spectra, scan_time_deltas, pixel_time_deltas)
        push!(processed_pixels, ProcessedPixel((x, y), mz, intensity))
    end

    return processed_pixels, (width, height)
end

function ExportImzml(target_file::String, pixels::Vector{ProcessedPixel}, dims::Tuple{Int, Int})
    @warn "ExportImzml is not yet implemented. The processed pixel data has been generated but not saved to .imzML/.ibd files."
    # TODO: Implement the logic to write the .imzML (XML) and .ibd (binary) files.
    # 1. Open target_file.imzML and target_file.ibd for writing.
    # 2. Write all m/z and intensity arrays sequentially to the .ibd file, tracking offsets and lengths.
    # 3. Write the .imzML XML structure, including:
    #    - Boilerplate headers
    #    - <fileDescription>, <referenceableParamGroupList>, <scanSettingsList>
    #    - A <spectrumList> with a <spectrum> for each pixel.
    #    - Each <spectrum> must contain cvParams for x/y coords and external data pointers to the .ibd file.
    return false
end

"""
    ImportMzmlFile(source_file::String, sync_file::String, target_file::String)

Main workflow function to convert a .mzML file to an .imzML file.

# Arguments
* `source_file`: Path to the input .mzML file.
* `sync_file`: Path to the synchronization text file.
* `target_file`: Path for the output .imzML file (the .ibd will be named accordingly).
"""
function ImportMzmlFile(source_file::String, sync_file::String, target_file::String)
    println("Step 1: Getting scan times from .mzML file...")
    scans = GetMzmlScanTime(source_file)

    println("Step 2: Matching acquisition times...")
    timing_matrix = MatchAcquireTime(sync_file, scans)

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
        println("Conversion failed. Exporting is not yet implemented.")
    end
    return success
end
