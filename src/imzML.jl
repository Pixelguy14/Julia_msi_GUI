
using Images, ImageFiltering, StatsBase, Statistics, CairoMakie, ColorSchemes, DataFrames, CSV, Printf, Dates
include("ParserHelpers.jl")

# --- Lazy Loading Data Structures ---

# Struct to hold metadata for a single spectrum, for lazy loading.
struct SpectrumMetadata
    x::Int32
    y::Int32
    mz_offset::Int64
    intensity_offset::Int64
    mz_count::Int32
    intensity_count::Int32
end

# Struct to hold the parsed imzML metadata and the handle to the .ibd file.
mutable struct ImzMLData
    ibd_handle::IO
    mz_format::Type
    intensity_format::Type
    spectra_metadata::Vector{SpectrumMetadata}
    width::Int
    height::Int
    mz_is_first::Bool # To know the order in the ibd file

    function ImzMLData(ibd_handle::IO, mz_format::Type, intensity_format::Type, spectra_metadata::Vector{SpectrumMetadata}, width::Int, height::Int, mz_is_first::Bool)
        obj = new(ibd_handle, mz_format, intensity_format, spectra_metadata, width, height, mz_is_first)
        finalizer(obj) do o
            if isopen(o.ibd_handle)
                close(o.ibd_handle)
            end
        end
        return obj
    end
end


# --- Extracted from imzML.jl ---

"""
This file provides a library for parsing `.imzML` and `.ibd` files in pure Julia.
It is intended to be included by a parent script.

Core Functions:
- `load_imzml`: The main function that orchestrates the parsing.
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
function axes_config_img(stream)
    tag = find_tag(stream, r"^\s*<(referenceableParamGroup )")
    value = get_attribute(tag.captures[1], "intensityArray")
    order = 1 + (value !== nothing)

    axis = Array{SpecDim,1}(undef, 2)
    axis[order] = configure_spec_dim(stream)

    find_tag(stream, r"^\s*<(referenceableParamGroup )")
    axis[xor(order, 3)] = configure_spec_dim(stream)
    return axis
end

"""
    get_img_dimensions(stream)

Reads the maximum X and Y dimensions and total spectrum count from the metadata.
"""
function get_img_dimensions(stream)
    find_tag(stream, r"^\s*<(scanSettings )")
    n = 2
    dim = [0, 0, 0]

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
function get_spectrum_tag_offset(stream)
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
function get_spectrum_attributes(stream, hIbd)
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


"""
    load_imzml(file_path::String)
Main function to parse an `.imzML`/.ibd file pair and prepare for lazy loading.
It parses the metadata from the `.imzML` file, including spectrum coordinates and
byte offsets for data in the `.ibd` file, but does not load the spectral data itself.

# Returns
- An `ImzMLData` struct containing the open `.ibd` file handle and all necessary metadata
  for on-demand data loading.
"""
function load_imzml(file_path::String)
    if !isfile(file_path)
        error("Provided path is not a file: $(file_path)")
    end

    ibd_path = replace(file_path, ".imzML" => ".ibd")
    if !isfile(ibd_path)
        error("Corresponding .ibd file not found for: $(file_path)")
    end

    stream = open(file_path)
    hIbd = open(ibd_path)

    # We don't use a try/finally block because we need to return the open hIbd handle.
    # The finalizer on the ImzMLData struct will be responsible for closing it.

    axis = axes_config_img(stream)
    imgDim = get_img_dimensions(stream)
    width, height, num_spectra = imgDim

    # Determine which format is for m/z and which for intensity from the referenceableParamGroups
    mz_config_idx = findfirst(a -> a.Axis == 1, axis)
    int_config_idx = findfirst(a -> a.Axis == 2, axis)
    mz_format = axis[mz_config_idx].Format
    intensity_format = axis[int_config_idx].Format

    # The original parsing logic is highly optimized and relies on fixed offsets within the XML structure.
    # We reuse it to get the initial .ibd file offset and the structural attributes of the <spectrum> tags.
    start_of_spectra_xml = position(stream)
    attr = get_spectrum_attributes(stream, hIbd)
    
    # get_spectrum_attributes moves the hIbd pointer to the start of the first spectrum's binary data.
    current_ibd_offset = position(hIbd)

    # We need to re-scan the spectra list from the beginning to get metadata for each spectrum.
    seek(stream, start_of_spectra_xml)

    spectra_metadata = Vector{SpectrumMetadata}(undef, num_spectra)
    
    # attr[3] is 4 if the first binary data array is intensity, 3 if it's m/z.
    mz_is_first = attr[3] == 3

    for k in 1:num_spectra
        # These skips are based on the structure of the first spectrum tag, assuming all are identical.
        skip(stream, attr[5]) # Skip to X coordinate value
        val_tag = find_tag(stream, r"value=\"(\d+)\"")
        x = parse(Int32, val_tag.captures[1])

        skip(stream, attr[6]) # Skip to Y coordinate value
        val_tag = find_tag(stream, r"value=\"(\d+)\"")
        y = parse(Int32, val_tag.captures[1])

        skip(stream, attr[7]) # Skip to array length value
        val_tag = find_tag(stream, r"value=\"(\d+)\"")
        nPoints = parse(Int32, val_tag.captures[1])

        mz_len_bytes = nPoints * sizeof(mz_format)
        int_len_bytes = nPoints * sizeof(intensity_format)

        local mz_offset, int_offset
        if mz_is_first
            mz_offset = current_ibd_offset
            int_offset = mz_offset + mz_len_bytes
        else
            int_offset = current_ibd_offset
            mz_offset = int_offset + int_len_bytes
        end

        spectra_metadata[k] = SpectrumMetadata(x, y, mz_offset, int_offset, nPoints, nPoints)
        
        # Advance the offset for the next spectrum's data.
        current_ibd_offset += mz_len_bytes + int_len_bytes
        
        skip(stream, attr[8]) # Skip to the end of the spectrum tag
    end
    
    # Close the .imzML file stream, but leave the .ibd stream open for lazy reading.
    close(stream)

    return ImzMLData(hIbd, mz_format, intensity_format, spectra_metadata, width, height, mz_is_first)
end

# --- End of content from imzML.jl ---


# --- Start of content from Imaging_Normalization.jl ---

# =============================================================================
#
# Image Slice Extraction
#
# =============================================================================

"""
    find_mass(mz_array, intensity_array, target_mass, tolerance)

Finds the intensity of the most intense peak within a mass tolerance window.
This modernized version is more robust than a simple binary search as it
correctly handles multiple peaks within the tolerance window.

# Returns
- The intensity (`Float64`) of the peak if found, otherwise `0.0`.
"""
function find_mass(mz_array, intensity_array, target_mass, tolerance)
    lower_bound = target_mass - tolerance
    upper_bound = target_mass + tolerance

    max_intensity = 0.0
    found = false

    # Iterate through the spectrum to find the highest intensity peak in the window
    for i in eachindex(mz_array)
        if lower_bound <= mz_array[i] <= upper_bound
            if intensity_array[i] > max_intensity
                max_intensity = intensity_array[i]
                found = true
            end
        end
    end

    return found ? max_intensity : 0.0
end

"""
    load_slices(folder, masses, tolerance)

Loads image slices for multiple masses from all `.imzML` files in a directory.
This function is optimized to read the binary data file (`.ibd`) only once per file.
It iterates through each spectrum, reads its data, and finds all target masses
before proceeding to the next spectrum, minimizing disk I/O.
"""
function load_slices(folder, masses, tolerance)
    files = filter(f -> endswith(f, ".imzML"), readdir(folder, join=true))
    if isempty(files)
        @warn "No .imzML files found in the specified directory: $folder"
        return (Array{Any}(undef, 0, 0), String[])
    end
    n_files = length(files)
    n_slices = length(masses)

    img_list = Array{Any}(undef, n_files, n_slices)
    names = String[]

    for (i, file) in enumerate(files)
        name = basename(file)
        push!(names, name)
        @info "Processing $(i)/$(n_files): $(name)"
        
        imzML_data = @time load_imzml(file) # This is fast, it's fine.

        # Create empty images for all slices for the current file
        current_file_slices = [zeros(Float64, imzML_data.width, imzML_data.height) for _ in 1:n_slices]

        hIbd = imzML_data.ibd_handle
        mz_format = imzML_data.mz_format
        int_format = imzML_data.intensity_format

        # Pre-allocate arrays to be reused for each spectrum
        mz_array = Array{mz_format}(undef, 0)
        intensity_array = Array{int_format}(undef, 0)

        # Iterate through each spectrum ONCE
        for meta in imzML_data.spectra_metadata
            # Resize arrays if the number of points in the spectrum has changed
            if length(mz_array) != meta.mz_count
                resize!(mz_array, meta.mz_count)
                resize!(intensity_array, meta.intensity_count)
            end

            # Read spectrum data ONCE from the .ibd file
            seek(hIbd, meta.mz_offset)
            read!(hIbd, mz_array)
            seek(hIbd, meta.intensity_offset)
            read!(hIbd, intensity_array)

            # Now, check for all masses of interest in this single spectrum
            for (j, mass) in enumerate(masses)
                intensity = find_mass(mz_array, intensity_array, mass, tolerance)
                if intensity > 0.0
                    if 1 <= meta.x <= imzML_data.width && 1 <= meta.y <= imzML_data.height
                        current_file_slices[j][meta.x, meta.y] = intensity
                    end
                end
            end
        end # end of spectra loop

        # Assign the generated images to the main list
        for j in 1:n_slices
            img_list[i, j] = current_file_slices[j]
        end

    end # end of files loop
    
    return (img_list, names)
end

# ============================================================================
#
#
# Image Processing and Normalization
#
# ============================================================================

function get_outlier_thres(img, prob=0.98)
    # DO NOT filter zeros. Use all pixel values like R does.
    int_values = vec(img)
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

function set_pixel_depth(img, bounds, depth)
    min_val, max_val = bounds
    bins = depth - 1
    
    if min_val >= max_val
        return zeros(UInt8, size(img))
    end
    
    # Create intensity bins
    range_vals = range(min_val, stop=max_val, length=depth)[2:depth]
    
    # Assign each pixel to a bin
    result = similar(img, UInt8) # Use similar to create an array of the same type and size
    for i in eachindex(img)
        if img[i] <= min_val
            result[i] = 0
        else
            bin_idx = findfirst(x -> img[i] <= x, range_vals)
            result[i] = bin_idx === nothing ? bins : bin_idx - 1
        end
    end
    
    return result
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

function median_filter(img)
    # 3x3 median filter implementation
    return mapwindow(median, img, (3, 3))
end

# ============================================================================
#
#
# Analysis and Visualization
#
# ============================================================================

"""
    display_statistics(slices)

Calculates and prints key statistics for each slice.
"""
function display_statistics(slices, names, masses)
    if isempty(slices)
        @warn "Cannot display statistics for empty slice list."
        return nothing
    end

    n_files, n_masses = size(slices)
    #stats_to_calc = Dict("Mean" => mean, "Max" => maximum, "Min" => minimum, "Sum" => sum, "Std" => std)
    stats_to_calc = Dict("Mean" => mean)
    all_dfs = Dict{String, DataFrame}()

    for (stat_name, stat_func) in stats_to_calc
        # Create a matrix to hold the statistic for each slice
        stat_matrix = zeros(Float64, n_files, n_masses)
        for i in 1:n_files, j in 1:n_masses
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

# =============================================================================
#
#
# Main Workflow
#
# ============================================================================

function load_config(path)
    config = Dict{String, Any}()
    if !isfile(path)
        @warn "Config file not found at $(path). Using default parameters."
        # Define defaults here in case the file is missing
        config["masses"] = [100]
        config["tolerance"] = 0.1
        config["color_depth"] = 255
        config["outlier_prob"] = 0.98
        return config
    end

    for line in eachline(path)
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue # Skip empty lines and comments
        end
        parts = split(line, '=', limit=2)
        if length(parts) != 2
            @warn "Skipping malformed line in config: $(line)"
            continue
        end
        key = strip(parts[1])
        value_str = strip(parts[2])

        try
            if key == "masses"
                values = [parse(Float64, s) for s in split(value_str, ',')]
                config[key] = values
            elseif key == "tolerance"
                config[key] = parse(Float64, value_str)
            elseif key == "color_depth"
                config[key] = parse(Int, value_str)
            elseif key == "outlier_prob"
                config[key] = parse(Float64, value_str)
            end
        catch e
            @error "Could not parse value for key '$(key)': $(value_str)"
        end
    end
    return config
end
