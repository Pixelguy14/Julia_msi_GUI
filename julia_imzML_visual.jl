# == Search functions ==
# Functions that recieve a list to update, and the current direction both as string for
# searching in the directory the position the list is going
function increment_image(current_image, image_list)
    if isempty(image_list)
        return nothing
    end
    current_index=findfirst(isequal(current_image), image_list)
    if current_index==nothing || current_index==length(image_list) || current_image ===""
        return image_list[length(image_list)]  # Return the current image if it's the last one or not found
    else
        return image_list[current_index + 1]  # Move to the next image
    end
end

function decrement_image(current_image, image_list)
    if isempty(image_list)
        return nothing
    end
    current_index=findfirst(isequal(current_image), image_list)
    if current_index==nothing || current_index==1 || current_image===""
        return image_list[1]  # Return the current image if it's the first one or not found
    else
        return image_list[current_index - 1]  # Move to the previous image
    end
end

## Plot Image functions
# Downsample an image matrix to a maximum dimension while preserving aspect ratio
function downsample_image(img_matrix, max_dim::Int)
    h, w = size(img_matrix)
    if h <= max_dim && w <= max_dim
        return img_matrix # No downsampling needed
    end

    aspect_ratio = w / h
    if w > h
        new_w = max_dim
        new_h = round(Int, max_dim / aspect_ratio)
    else
        new_h = max_dim
        new_w = round(Int, max_dim * aspect_ratio)
    end

    # imresize from Images.jl is perfect for this
    return imresize(img_matrix, (new_h, new_w))
end

# loadImgPlot recieves the local directory of the image as a string,
# returns the layout and data for the heatmap plotly plot
# this function loads the image into a plot
function loadImgPlot(interfaceImg::String)
    # Load the image
    cleaned_img=replace(interfaceImg, r"\?.*" => "")
    cleaned_img=lstrip(cleaned_img, '/')
    var=joinpath("./public", cleaned_img)
    img=load(var)
    # Convert to grayscale
    img_gray=Gray.(img)
    img_array=Array(img_gray)
    elevation=Float32.(Array(img_array)) ./ 255.0
    # Get the X, Y coordinates of the image
    height, width=size(img_array)
    X=collect(1:width)
    Y=collect(1:height)

    # Create the layout
    layout=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=14,
                color="black"
            )
        ),
        xaxis=PlotlyBase.attr(
            visible=false,
            scaleanchor="y",
            range=[0, width]
        ),
        yaxis=PlotlyBase.attr(
            visible=false,
            range=[-height, 0]
        ),
        margin=attr(l=0,r=0,t=0,b=0,pad=0)
    )

    # Create the trace for the image
    trace=PlotlyBase.heatmap(
        z=elevation,
        x=X,
        y=-Y,
        name="",
        hoverinfo="x+y",
        showlegend=false,
        colorscale="Viridis",
        showscale=false,
        colorbar=attr(
            title=attr(
                text="Intensity",
                font=attr(
                    size=14,
                    color="black"
                ),
                side="right"
            ),
            ticks="outside",
            ticklen=2,
            tickwidth=0.5,
            nticks=5,
            tickformat=".2g"
        )
    )

    plotdata=[trace]
    plotlayout=layout
    return plotdata, plotlayout, width, height
end

# loadImgPlot recieves the local directory of the image as a string, the local directory o the overlay image
# and the transparency its required to have. Returns the layout and data for the heatmap plotly plot
# this function loads the image into a plot
function loadImgPlot(interfaceImg::String, overlayImg::String, imgTrans::Float64)
    timestamp=string(time_ns())
    # Load the main image
    cleaned_img = replace(interfaceImg, r"\?.*" => "")
    cleaned_img = lstrip(cleaned_img, '/')
    var = joinpath("./public", cleaned_img)
    img = load(var)
    # Convert to grayscale
    img_gray = Gray.(img)
    img_array = Array(img_gray)
    elevation = Float32.(Array(img_array)) ./ 255.0
    # Get the X, Y coordinates of the image
    height, width = size(img_array)
    X = collect(1:width)
    Y = collect(1:height)

    # Create the layout with overlay image
    layoutImg = PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=14,
                color="black"
            )
        ),
        images = [attr(
            source = "$(overlayImg)?t=$(timestamp)",
            xref = "x",
            yref = "y",
            x = 0,
            y = 0,
            sizex = width,
            sizey = -height,
            sizing = "stretch",
            opacity = imgTrans,
            layer = "above"  # Place the overlay image in the foreground
        )],
        xaxis = PlotlyBase.attr(
            visible = false,
            scaleanchor = "y",
            range = [0, width]
        ),
        yaxis = PlotlyBase.attr(
            visible = false,
            range = [-height,0]
        ),
        margin = attr(l = 0, r = 0, t = 0, b = 0, pad = 0)
    )

    # Create the trace for the main image
    trace = PlotlyBase.heatmap(
        z = elevation,
        x = X,
        y = -Y,
        name = "",
        hoverinfo = "x+y",
        showlegend = false,
        colorscale = "Viridis",
        showscale = false
    )

    plotdata = [trace]
    plotlayout = layoutImg
    return plotdata, plotlayout, width, height
end

# loadContourPlot recieves the local directory of the image as a string,
# returns the layout and data for the contour plotly plot
# this function loads the image and applies a gaussian filter 
# to smoothen it and loads it into a plot
function loadContourPlot(interfaceImg::String)
    # Load the image
    cleaned_img=replace(interfaceImg, r"\?.*" => "")
    cleaned_img=lstrip(cleaned_img, '/')
    var=joinpath("./public", cleaned_img)
    img=load(var)
    img_gray=Gray.(img)
    img_array=Array(img_gray)
    elevation=Float32.(Array(img_array))./ 255.0 # Normalize between 0 and 1
    
    # Smooth the image
    sigma=3.0
    kernel=Kernel.gaussian(sigma)
    elevation_smoothed=imfilter(elevation, kernel)

    # --- DOWNSAMPLING FOR PERFORMANCE ---
    elevation_smoothed = downsample_image(elevation_smoothed, 512)
    # ---
    
    # Create the X, Y meshgrid coordinates
    x=1:size(elevation_smoothed, 2)
    y=1:size(elevation_smoothed, 1)
    X=repeat(reshape(x, 1, length(x)), length(y), 1)
    Y=repeat(reshape(y, length(y), 1), 1, length(x))

    # Define tick values and text for colorbars
    min_val = minimum(elevation_smoothed)
    max_val = maximum(elevation_smoothed)
    tickV = range(min_val, stop=max_val, length=8)
    tickT = log_tick_formatter(collect(tickV))
                
    layout=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="2D topographic map of $cleaned_img (downsampled)",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        xaxis=PlotlyBase.attr(
            visible=false,
            scaleanchor="y"
        ),
        yaxis=PlotlyBase.attr(
            visible=false
        ),
        margin=attr(l=0,r=0,t=100,b=0,pad=0)
    )
    trace=PlotlyBase.contour(
        z=elevation_smoothed,
        x=X[1, :],  # Use the first row
        y=-Y[:, 1],  # Use the first column
        contours_coloring="Viridis", 
        colorscale="Viridis",
        colorbar = attr(
            tickvals = tickV,
            ticktext = tickT,
            tickmode = "array"
        )
    )
    plotdata=[trace]
    plotlayout=layout
    return plotdata, plotlayout
end

# loadSurfacePlot recieves the local directory of the image as a string,
# returns the layout and data for the surface plotly plot
# this function loads the image and applies a gaussian filter 
# to smoothen it and loads it into a 3D plot
function loadSurfacePlot(interfaceImg::String)
    # Load the image
    cleaned_img=replace(interfaceImg, r"\?.*" => "")
    cleaned_img=lstrip(cleaned_img, '/')
    var=joinpath("./public", cleaned_img)
    img=load(var)
    img_gray=Gray.(img) # Convert to grayscale
    img_array=Array(img_gray)
    elevation=Float32.(Array(img_array)) ./ 255.0 # Normalize between 0 and 1
    
    # Smooth the image 
    sigma=3.0
    kernel=Kernel.gaussian(sigma)
    elevation_smoothed=imfilter(elevation, kernel)

    # --- DOWNSAMPLING FOR PERFORMANCE ---
    elevation_smoothed = downsample_image(elevation_smoothed, 256)
    # ---
    
    # Create the X, Y meshgrid coordinates 
    x=1:size(elevation_smoothed, 2) 
    y=1:size(elevation_smoothed, 1)
    X=repeat(reshape(x, 1, length(x)), length(y), 1)
    Y=repeat(reshape(y, length(y), 1), 1, length(x))

    # Define tick values and text for colorbars
    min_val = minimum(elevation_smoothed)
    max_val = maximum(elevation_smoothed)
    tickV = range(min_val, stop=max_val, length=8)
    tickT = log_tick_formatter(collect(tickV))

    # Calculate the number of ticks and aspect ratio for the 3d plot
    x_nticks=min(20, length(x))
    y_nticks=min(20, length(y))
    z_nticks=5
    aspect_ratio=attr(x=1, y=length(y) / length(x), z=0.5)
    # Define the layout for the 3D plot
    layout3D=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="3D surface plot of $cleaned_img (downsampled)",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        scene=attr(
            xaxis_nticks=x_nticks,
            yaxis_nticks=y_nticks,
            zaxis_nticks=z_nticks, 
            camera=attr(eye=attr(x=0, y=1, z=0.5)), 
            aspectratio=aspect_ratio
        ),
        margin=attr(l=0,r=0,t=120,b=0,pad=0)
    )
    # Transpose the elevation_smoothed array if Y axis is longer than X axis to fix chopping
    elevation_smoothed=transpose(elevation_smoothed)
    if size(elevation_smoothed, 1) < size(elevation_smoothed, 2)
        Y=-Y
    else
        X=-X
    end

    trace3D=PlotlyBase.surface(
        x=X[1, :], 
        y=Y[:, 1], 
        z=elevation_smoothed,
        contours_z=attr(
            show=true,
            usecolormap=true,
            highlightcolor="limegreen",
            project_z=true
        ), 
        colorscale="Viridis",
        colorbar = attr(
            tickvals = tickV,
            ticktext = tickT,
            nticks=8
        )
    )
    plotdata=[trace3D]
    plotlayout=layout3D
    return plotdata, plotlayout
end

# This function recieves the x and y coords currently selected, and the dimentions of
# the image to create two traces that will display in a cross section
function crossLinesPlot(x, y, maxwidth, maxheight)
    # Define the coordinates for the two lines
    l1_x=[0, maxwidth]
    l1_y=[y, y]
    l2_x=[x, x]
    l2_y=[0, maxheight]

    # Create the line traces
    trace1=PlotlyBase.scatter(x=l1_x, y=l1_y, mode="lines",line=attr(color="red", width=0.5),name="Line X",showlegend=false)
    trace2=PlotlyBase.scatter(x=l2_x, y=l2_y, mode="lines",line=attr(color="red", width=0.5),name="Line Y",showlegend=false)

    return trace1, trace2
end

# This function is used for giving colorbar values a visual format
# that shortens long values giving them scientific notation
function log_tick_formatter(values::Vector{Float64})
    # Initialize exponents dictionary
    exponents=zeros(Int, length(values))
    formValues=zeros(Float64, length(values))
    for i in 1:length(values)
        value = values[i]
        if value >= 1000  # positive formatting for notation
            while value >= 1000
                value /= 10
                exponents[i] += 1
            end
        elseif value > 0 && value < 1  # negative formatting for notation
            while value < 1
                value *= 10
                exponents[i] -= 1
            end
        end
        formValues[i]=value
    end
    return map((v, e) -> e == 0 ? "$(round(v, sigdigits=2))" : "$(round(v, sigdigits=2))x10" * Makie.UnicodeFun.to_superscript(e), formValues, exponents)
end

function generate_colorbar_image(slice_data::AbstractMatrix, color_levels::Int, output_path::String; use_triq::Bool=false, triq_prob::Float64=0.98)
    # 1. Determine bounds based on whether TrIQ is used
    min_val, max_val = if use_triq
        MSI_src.get_outlier_thres(slice_data, triq_prob)
    else
        extrema(slice_data)
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

# meanSpectrumPlot recieves the local directory of the image as a string,
# returns the layout and data for the surface plotly plot
# this function loads the spectra data and makes a mean to display
# its values in the spectrum plot
function meanSpectrumPlot(data::MSIData, dataset_name::String="")
    title_text = isempty(dataset_name) ? "Average Spectrum Plot" : "Average Spectrum for: $dataset_name"
    layout = PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text=title_text,
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        hovermode="closest",
        xaxis=PlotlyBase.attr(
            title="<i>m/z</i>",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Average Intensity",
            showgrid=true,
            tickformat=".3g"
        ),
        margin=attr(l=0, r=0, t=120, b=0, pad=0)
    )

    # Use the new, efficient function from the backend
    xSpectraMz, ySpectraMz = get_average_spectrum(data)

    if isempty(xSpectraMz)
        @warn "Average spectrum is empty."
        trace = PlotlyBase.scatter(x=Float64[], y=Float64[])
    else
        trace = PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, marker=attr(size=1, color="blue", opacity=0.5), name="Average", hoverinfo="x",hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    end


    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end

function xySpectrumPlot(data::MSIData, xCoord::Int, yCoord::Int, imgWidth::Int, imgHeight::Int, dataset_name::String="")
    local mz::AbstractVector, intensity::AbstractVector
    local plot_title::String

    is_imaging = data.source isa ImzMLSource

    if is_imaging
        # For imaging data, use (X, Y) coordinates
        x = clamp(xCoord, 1, imgWidth)
        y = clamp(yCoord, 1, imgHeight)
        
        process_spectrum(data, Int(x), Int(y)) do recieved_mz, recieved_intensity
            mz = recieved_mz
            intensity = recieved_intensity
        end
        base_title = "Spectrum at ($x, $y)"
    else
        # For non-imaging data, treat xCoord as the spectrum index
        index = clamp(xCoord, 1, length(data.spectra_metadata))
        
        process_spectrum(data, index)  do recieved_mz, recieved_intensity
            mz = recieved_mz
            intensity = recieved_intensity
        end
        base_title = "Spectrum #$index"
    end

    plot_title = isempty(dataset_name) ? base_title : "$base_title for: $dataset_name"

    layout = PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text=plot_title,
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        hovermode="closest",
        xaxis=PlotlyBase.attr(
            title="<i>m/z</i>",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Intensity",
            showgrid=true,
            tickformat=".3g"
        ),
        margin=attr(l=0, r=0, t=120, b=0, pad=0)
    )
    
    # Downsample for plotting performance
    mz_down, int_down = MSI_src.downsample_spectrum(mz, intensity)

    trace = PlotlyBase.scatter(x=mz_down, y=int_down, marker=attr(size=1, color="blue", opacity=0.5), name="Spectrum", hoverinfo="x", hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    
    plotdata = [trace]
    plotlayout = layout
    
    return plotdata, plotlayout, mz, intensity
end

function sumSpectrumPlot(data::MSIData, dataset_name::String="")
    title_text = isempty(dataset_name) ? "Total Spectrum Plot" : "Total Spectrum for: $dataset_name"
    layout = PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text=title_text,
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        hovermode="closest",
        xaxis=PlotlyBase.attr(
            title="<i>m/z</i>",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Total Intensity",
            showgrid=true,
            tickformat=".3g"
        ),
        margin=attr(l=0, r=0, t=120, b=0, pad=0)
    )

    # Use the get_total_spectrum function from the backend
    xSpectraMz, ySpectraMz = get_total_spectrum(data)

    if isempty(xSpectraMz)
        @warn "Total spectrum is empty."
        trace = PlotlyBase.scatter(x=Float64[], y=Float64[])
    else
        trace = PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, marker=attr(size=1, color="blue", opacity=0.5), name="Total", hoverinfo="x",hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    end

    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end

function warmup_init()
    @async begin
        println("Pre-compiling functions at startup...")
        
        # Create a dummy MSIData object to be used for pre-compilation
        # dummy_source = ImzMLSource("dummy.ibd", Float32, Float32)
        # dummy_meta = MSI_src.SpectrumMetadata(0,0,"",MSI_src.UNKNOWN, MSI_src.SpectrumAsset(Float32,false,0,0,:mz), MSI_src.SpectrumAsset(Float32,false,0,0,:intensity))
        # dummy_msi_data = MSIData(dummy_source, [dummy_meta], (1,1), zeros(Int,1,1), 0)

        # Pre-compile functions from btnSearch
        # try OpenMSIData("dummy.imzML") catch end
        # try precompute_analytics(dummy_msi_data) catch end

        # Pre-compile functions from mainProcess
        # try get_mz_slice(dummy_msi_data, 1.0, 1.0) catch end
        try TrIQ(zeros(10,10), 256, 0.98) catch end
        try quantize_intensity(zeros(10,10), 256) catch end
        
        dummy_bmp_path = joinpath("public", "dummy.bmp")
        dummy_png_path = joinpath("public", "dummy.png")
        try 
            save_bitmap(dummy_bmp_path, zeros(UInt8, 10, 10), ViridisPalette)
            loadImgPlot("/dummy.bmp")
            generate_colorbar_image(zeros(10,10), 256, dummy_png_path)
        catch e
            @warn "Pre-compilation step failed (this is expected if dummy files can't be created/read)"
        finally
            rm(dummy_bmp_path, force=true)
            rm(dummy_png_path, force=true)
        end

        println("Pre-compilation finished.")
    end
end

function load_registry(registry_path)
    if isfile(registry_path)
        try
            return JSON.parsefile(registry_path, dicttype=Dict{String, Any})
        catch e
            @error "Failed to parse registry.json: $e"
            return Dict{String, Any}()
        end
    end
    return Dict{String, Any}()
end

function extract_metadata(msi_data::MSIData, source_path::String)
    df = msi_data.spectrum_stats_df
    if df === nothing
        # This can happen if precompute_analytics hasn't been run
        # We can still return basic info
        return Dict(
            "summary" => [
                Dict("parameter" => "File Name", "value" => basename(source_path)),
                Dict("parameter" => "Number of Spectra", "value" => length(msi_data.spectra_metadata)),
                Dict("parameter" => "Image Dimensions", "value" => "$(msi_data.image_dims[1]) x $(msi_data.image_dims[2])"),
            ],
            "global_min_mz" => nothing,
            "global_max_mz" => nothing
        )
    end

    summary_stats = [
        Dict("parameter" => "File Name", "value" => basename(source_path)),
        Dict("parameter" => "Number of Spectra", "value" => length(msi_data.spectra_metadata)),
        Dict("parameter" => "Image Dimensions", "value" => "$(msi_data.image_dims[1]) x $(msi_data.image_dims[2])"),
        Dict("parameter" => "Global Min m/z", "value" => @sprintf("%.4f", msi_data.global_min_mz)),
        Dict("parameter" => "Global Max m/z", "value" => @sprintf("%.4f", msi_data.global_max_mz)),
        Dict("parameter" => "Mean TIC", "value" => @sprintf("%.2e", mean(df.TIC))),
        Dict("parameter" => "Mean BPI", "value" => @sprintf("%.2e", mean(df.BPI))),
        Dict("parameter" => "Mean # Points", "value" => @sprintf("%.1f", mean(df.NumPoints))),
    ]

    if hasproperty(df, :Mode)
        centroid_count = count(==(MSI_src.CENTROID), df.Mode)
        profile_count = count(==(MSI_src.PROFILE), df.Mode)
        unknown_count = count(==(MSI_src.UNKNOWN), df.Mode)
        
        push!(summary_stats, Dict("parameter" => "Centroid Spectra", "value" => string(centroid_count)))
        push!(summary_stats, Dict("parameter" => "Profile Spectra", "value" => string(profile_count)))
        if unknown_count > 0
            push!(summary_stats, Dict("parameter" => "Unknown Mode Spectra", "value" => string(unknown_count)))
        end
    end

    return Dict(
        "summary" => summary_stats,
        "global_min_mz" => msi_data.global_min_mz,
        "global_max_mz" => msi_data.global_max_mz
    )
end

function update_registry(registry_path, dataset_name, source_path, metadata=nothing, is_imzML=false)
    registry = load_registry(registry_path)
    entry = Dict{String, Any}( # Explicitly type the dictionary to allow mixed value types
        "source_path" => source_path, 
        "processed_date" => string(now()),
        "is_imzML" => is_imzML
    )
    if metadata !== nothing
        entry["metadata"] = metadata
    end
    registry[dataset_name] = entry
    
    try
        open(registry_path, "w") do f
            JSON.print(f, registry, 4)
        end
    catch e
        @error "Failed to write to registry.json: $e"
    end
end

function process_file_safely(file_path, masses, params, progress_message_ref, overall_progress_ref)
    local_msi_data = nothing
    dataset_name = replace(basename(file_path), r"\.imzML$"i => "")
    output_dir = joinpath("public", dataset_name)
    println("Processing: $dataset_name -> $output_dir")
    
    try
        # --- Load Data ---
        progress_message_ref = "Loading: $(basename(file_path))"
        local_msi_data = OpenMSIData(file_path)
        if !(local_msi_data.source isa ImzMLSource)
            @warn "Skipping non-imzML file: $(basename(file_path))"
            return (false, "Skipped: Not an imzML file")
        end

        # --- Generate Slices (this will call precompute_analytics if needed) ---
        progress_message_ref = "Generating $(length(masses)) slices for $(dataset_name)..."
        slice_dict = get_multiple_mz_slices(local_msi_data, masses, params.tolerance)

        # --- Extract metadata *after* it has been computed ---
        metadata = extract_metadata(local_msi_data, file_path)

        # --- Save Slices ---
        mkpath(output_dir)  # Ensure output directory exists
        
        for (mass_idx, mass) in enumerate(masses)
            progress_message_ref = "File $(params.fileIdx)/$(params.nFiles): Saving slice for m/z=$mass"
            
            slice = slice_dict[mass]
            text_nmass = replace(string(mass), "." => "_")
            bitmap_filename = params.triqE ? "TrIQ_$(text_nmass).bmp" : "MSI_$(text_nmass).bmp"
            colorbar_filename = params.triqE ? "colorbar_TrIQ_$(text_nmass).png" : "colorbar_MSI_$(text_nmass).png"

            if all(iszero, slice)
                sliceQuant = zeros(UInt8, size(slice))
                @warn "No intensity data for m/z = $mass in $(dataset_name)"
            else
                sliceQuant = params.triqE ? TrIQ(slice, params.colorL, params.triqP) : quantize_intensity(slice, params.colorL)
                if params.medianF
                    sliceQuant = round.(UInt8, median_filter(sliceQuant))
                end
            end

            save_bitmap(joinpath(output_dir, bitmap_filename), sliceQuant, ViridisPalette)
            if !all(iszero, slice)
                generate_colorbar_image(slice, params.colorL, joinpath(output_dir, colorbar_filename); use_triq=params.triqE, triq_prob=params.triqP)
            end
        end
        
        is_imzML = local_msi_data.source isa ImzMLSource
        update_registry(params.registry, dataset_name, file_path, metadata, is_imzML)
        return (true, "")

    catch e
        @error "File processing failed" file=file_path exception=(e, catch_backtrace())
        return (false, "File: $(basename(file_path)) - $(sprint(showerror, e))")
    finally
        if local_msi_data !== nothing
            # Cleanup
        end
        local_msi_data = nothing
        GC.gc(true)
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0)
        end
    end
end