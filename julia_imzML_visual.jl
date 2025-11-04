# julia_imzML_visual.jl

"""
    increment_image(current_image, image_list)

Finds the next image in a list.

# Arguments
- `current_image`: The current image file name.
- `image_list`: The list of available image file names.

# Returns
- The file name of the next image, or the last image if the current one is the last or not found.
- `nothing` if `image_list` is empty.
"""
function increment_image(current_image, image_list)
    if isempty(image_list)
        return nothing
    end
    current_index = findfirst(isequal(current_image), image_list)
    if current_index === nothing || current_index == length(image_list) || current_image == ""
        return image_list[end]  # Return the last image if current is not found or is the last
    else
        return image_list[current_index + 1]  # Move to the next image
    end
end

"""
    decrement_image(current_image, image_list)

Finds the previous image in a list.

# Arguments
- `current_image`: The current image file name.
- `image_list`: The list of available image file names.

# Returns
- The file name of the previous image, or the first image if the current one is the first or not found.
- `nothing` if `image_list` is empty.
"""
function decrement_image(current_image, image_list)
    if isempty(image_list)
        return nothing
    end
    current_index = findfirst(isequal(current_image), image_list)
    if current_index === nothing || current_index == 1 || current_image == ""
        return image_list[1]  # Return the first image if current is not found or is the first
    else
        return image_list[current_index - 1]  # Move to the previous image
    end
end

"""
    downsample_image(img_matrix, max_dim::Int)

Downsamples an image matrix to a maximum dimension while preserving aspect ratio.

# Arguments
- `img_matrix`: The image matrix to downsample.
- `max_dim`: The maximum dimension (width or height) for the downsampled image.

# Returns
- The downsampled image matrix.
"""
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

"""
    loadImgPlot(interfaceImg::String)

Loads an image and creates a Plotly heatmap.

# Arguments
- `interfaceImg`: The path to the image file relative to the "public" directory.

# Returns
- `plotdata`: A vector containing the Plotly trace.
- `plotlayout`: The Plotly layout for the plot.
- `width`: The width of the loaded image.
- `height`: The height of the loaded image.
"""
function loadImgPlot(interfaceImg::String)
    # Load the image
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
    X = 1:width
    Y = 1:height

    # Create the layout
    layout = PlotlyBase.Layout(
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
        margin=attr(l=0, r=0, t=0, b=0, pad=0)
    )

    # Create the trace for the image
    trace = PlotlyBase.heatmap(
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

    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, width, height
end

"""
    loadImgPlot(interfaceImg::String, overlayImg::String, imgTrans::Float64)

Loads a main image and overlays a second image on top, creating a Plotly heatmap.

# Arguments
- `interfaceImg`: Path to the main image file.
- `overlayImg`: Path to the overlay image file.
- `imgTrans`: Transparency level for the overlay image (0.0 to 1.0).

# Returns
- `plotdata`: A vector containing the Plotly trace for the main image.
- `plotlayout`: The Plotly layout, including the overlay image.
- `width`: The width of the main image.
- `height`: The height of the main image.
"""
function loadImgPlot(interfaceImg::String, overlayImg::String, imgTrans::Float64)
    timestamp = string(time_ns())
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
    X = 1:width
    Y = 1:height

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
        images=[attr(
            source="$(overlayImg)?t=$(timestamp)",
            xref="x",
            yref="y",
            x=0,
            y=0,
            sizex=width,
            sizey=-height,
            sizing="stretch",
            opacity=imgTrans,
            layer="above"  # Place the overlay image in the foreground
        )],
        xaxis=PlotlyBase.attr(
            visible=false,
            scaleanchor="y",
            range=[0, width]
        ),
        yaxis=PlotlyBase.attr(
            visible=false,
            range=[-height, 0]
        ),
        margin=attr(l=0, r=0, t=0, b=0, pad=0)
    )

    # Create the trace for the main image
    trace = PlotlyBase.heatmap(
        z=elevation,
        x=X,
        y=-Y,
        name="",
        hoverinfo="x+y",
        showlegend=false,
        colorscale="Viridis",
        showscale=false
    )

    plotdata = [trace]
    plotlayout = layoutImg
    return plotdata, plotlayout, width, height
end

"""
    loadContourPlot(interfaceImg::String)

Loads an image, smooths it, and creates a Plotly contour plot.

# Arguments
- `interfaceImg`: Path to the image file.

# Returns
- `plotdata`: A vector containing the Plotly contour trace.
- `plotlayout`: The Plotly layout for the plot.
"""
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

"""
    loadSurfacePlot(interfaceImg::String)

Loads an image, smooths it, and creates a 3D Plotly surface plot.

# Arguments
- `interfaceImg`: Path to the image file.

# Returns
- `plotdata`: A vector containing the Plotly surface trace.
- `plotlayout`: The Plotly layout for the 3D plot.
"""
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

"""
    crossLinesPlot(x, y, maxwidth, maxheight)

Creates two line traces for a crosshair indicator on a plot.

# Arguments
- `x`: The x-coordinate of the crosshair center.
- `y`: The y-coordinate of the crosshair center.
- `maxwidth`: The width of the plot area.
- `maxheight`: The height of the plot area.

# Returns
- `trace1`: The horizontal line trace.
- `trace2`: The vertical line trace.
"""
function crossLinesPlot(x, y, maxwidth, maxheight)
    # Define the coordinates for the two lines
    l1_x = [0, maxwidth]
    l1_y = [y, y]
    l2_x = [x, x]
    l2_y = [0, maxheight]

    # Create the line traces
    trace1 = PlotlyBase.scatter(x=l1_x, y=l1_y, mode="lines", line=attr(color="red", width=0.5), name="Line X", showlegend=false)
    trace2 = PlotlyBase.scatter(x=l2_x, y=l2_y, mode="lines", line=attr(color="red", width=0.5), name="Line Y", showlegend=false)

    return trace1, trace2
end

"""
    log_tick_formatter(values::Vector{Float64})

Formats a vector of numbers into strings with a custom scientific notation for use as tick labels.
For example, 1000 becomes "100x10¹" and 0.01 becomes "1.0x10⁻²".

# Arguments
- `values`: A vector of `Float64` values to format.

# Returns
- A vector of formatted strings.
"""
function log_tick_formatter(values::Vector{Float64})
    # Initialize exponents dictionary
    exponents = zeros(Int, length(values))
    formValues = zeros(Float64, length(values))
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
        formValues[i] = value
    end
    return map((v, e) -> e == 0 ? "$(round(v, sigdigits=2))" : "$(round(v, sigdigits=2))x10" * Makie.UnicodeFun.to_superscript(e), formValues, exponents)
end

"""
    meanSpectrumPlot(data::MSIData, dataset_name::String="")

Generates a plot of the mean spectrum from MSIData.

# Arguments
- `data`: The `MSIData` object.
- `dataset_name`: Optional name of the dataset for the plot title.

# Returns
- `plotdata`: A vector containing the Plotly trace.
- `plotlayout`: The Plotly layout for the plot.
- `xSpectraMz`: The m/z values of the spectrum.
- `ySpectraMz`: The intensity values of the spectrum.
"""
function meanSpectrumPlot(data::MSIData, dataset_name::String=""; mask_path::Union{String, Nothing}=nothing)
    # Determine base title based on mask usage
    base_title = if mask_path !== nothing
        "Masked Average Spectrum"
    else
        "Average Spectrum"
    end
    
    title_text = isempty(dataset_name) ? base_title : "$base_title for: $dataset_name"
    
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
    xSpectraMz, ySpectraMz = get_average_spectrum(data, mask_path=mask_path)

    if isempty(xSpectraMz) || isempty(ySpectraMz)
        @warn "Average spectrum is empty."
        trace = PlotlyBase.stem(x=Float64[], y=Float64[])
        # Update title to indicate empty spectrum
        layout.title.text = "Empty " * layout.title.text
    else
        trace = PlotlyBase.stem(x=xSpectraMz, y=ySpectraMz, marker=attr(size=1, color="blue", opacity=0.5), name="Average", hoverinfo="x", hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    end

    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end

"""
    xySpectrumPlot(data::MSIData, xCoord::Int, yCoord::Int, imgWidth::Int, imgHeight::Int, dataset_name::String="")

Generates a plot for the spectrum at a specific coordinate (for imaging data) or index (for non-imaging data).

# Arguments
- `data`: The `MSIData` object.
- `xCoord`: The x-coordinate or spectrum index.
- `yCoord`: The y-coordinate (used for imaging data).
- `imgWidth`: The width of the MSI image.
- `imgHeight`: The height of the MSI image.
- `dataset_name`: Optional name of the dataset for the plot title.

# Returns
- `plotdata`: A vector containing the Plotly trace.
- `plotlayout`: The Plotly layout for the plot.
- `mz`: The m/z values of the spectrum.
- `intensity`: The intensity values of the spectrum.
"""
function xySpectrumPlot(data::MSIData, xCoord::Int, yCoord::Int, imgWidth::Int, imgHeight::Int, dataset_name::String=""; mask_path::Union{String, Nothing}=nothing)
    local mz::AbstractVector, intensity::AbstractVector
    local plot_title::String

    is_imaging = data.source isa ImzMLSource

    if is_imaging
        x = clamp(xCoord, 1, imgWidth)
        y = clamp(yCoord, 1, imgHeight)

        if mask_path !== nothing
            mask_matrix = load_and_prepare_mask(mask_path, (imgWidth, imgHeight))
            if !mask_matrix[y, x]
                @warn "Coordinate ($x, $y) is outside the specified mask."
                mz, intensity = Float64[], Float64[]
                base_title = "Spectrum Outside Mask at ($x, $y)"
            else
                process_spectrum(data, Int(x), Int(y)) do recieved_mz, recieved_intensity
                    mz = recieved_mz
                    intensity = recieved_intensity
                end
                base_title = "Masked Spectrum at ($x, $y)"
            end
        else
            process_spectrum(data, Int(x), Int(y)) do recieved_mz, recieved_intensity
                mz = recieved_mz
                intensity = recieved_intensity
            end
            base_title = "Spectrum at ($x, $y)"
        end
    else
        # For non-imaging data, treat xCoord as the spectrum index
        index = clamp(xCoord, 1, length(data.spectra_metadata))

        process_spectrum(data, index) do recieved_mz, recieved_intensity
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

    trace = PlotlyBase.stem(x=mz_down, y=int_down, marker=attr(size=1, color="blue", opacity=0.5), name="Spectrum", hoverinfo="x", hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")

    plotdata = [trace]
    plotlayout = layout

    return plotdata, plotlayout, mz, intensity
end

"""
    sumSpectrumPlot(data::MSIData, dataset_name::String="")

Generates a plot of the total (summed) spectrum from MSIData.

# Arguments
- `data`: The `MSIData` object.
- `dataset_name`: Optional name of the dataset for the plot title.

# Returns
- `plotdata`: A vector containing the Plotly trace.
- `plotlayout`: The Plotly layout for the plot.
- `xSpectraMz`: The m/z values of the spectrum.
- `ySpectraMz`: The intensity values of the spectrum.
"""
function sumSpectrumPlot(data::MSIData, dataset_name::String=""; mask_path::Union{String, Nothing}=nothing)
    # Determine base title based on mask usage
    base_title = if mask_path !== nothing
        "Masked Total Spectrum"
    else
        "Total Spectrum"
    end
    
    title_text = isempty(dataset_name) ? base_title : "$base_title for: $dataset_name"
    
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
    xSpectraMz, ySpectraMz, num_spectra = get_total_spectrum(data, mask_path=mask_path)

    if isempty(xSpectraMz) || isempty(ySpectraMz)
        @warn "Total spectrum is empty."
        trace = PlotlyBase.stem(x=Float64[], y=Float64[])
        # Update title to indicate empty spectrum
        layout.title.text = "Empty " * layout.title.text
    else
        trace = PlotlyBase.stem(x=xSpectraMz, y=ySpectraMz, marker=attr(size=1, color="blue", opacity=0.5), name="Total", hoverinfo="x", hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    end

    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end

"""
    warmup_init()

Performs pre-compilation of key functions at application startup to reduce first-use latency.
This is run asynchronously and should not block application startup.
"""
function warmup_init()
    @async begin
        println("Pre-compiling functions at startup...")

        # Pre-compile image processing and plotting functions
        try
            TrIQ(zeros(10, 10), 256, 0.98)
        catch
        end
        try
            quantize_intensity(zeros(10, 10), 256)
        catch
        end

        dummy_bmp_path = joinpath("public", "dummy.bmp")
        dummy_png_path = joinpath("public", "dummy.png")
        try
            save_bitmap(dummy_bmp_path, zeros(UInt8, 10, 10), ViridisPalette)
            loadImgPlot("/dummy.bmp")
            generate_colorbar_image(zeros(10, 10), 256, dummy_png_path)
        catch e
            @warn "Pre-compilation step failed (this is expected if dummy files can't be created/read)"
        finally
            rm(dummy_bmp_path, force=true)
            rm(dummy_png_path, force=true)
        end

        println("Pre-compilation finished.")
    end
end

"""
    load_registry(registry_path)

Loads the dataset registry from a JSON file.

# Arguments
- `registry_path`: Path to the `registry.json` file.

# Returns
- A dictionary containing the registry data. Returns an empty dictionary if the file doesn't exist or fails to parse.
"""
function load_registry(registry_path)
    if isfile(registry_path)
        try
            return JSON.parsefile(registry_path, dicttype=Dict{String,Any})
        catch e
            @error "Failed to parse registry.json: $e"
            return Dict{String,Any}()
        end
    end
    return Dict{String,Any}()
end

"""
    extract_metadata(msi_data::MSIData, source_path::String)

Extracts key metadata from an `MSIData` object for display.

# Arguments
- `msi_data`: The `MSIData` object.
- `source_path`: The path to the source data file.

# Returns
- A dictionary containing summary statistics and metadata.
"""
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

"""
    update_registry(registry_path, dataset_name, source_path, metadata=nothing, is_imzML=false)

Adds or updates an entry in the dataset registry JSON file.

# Arguments
- `registry_path`: Path to the `registry.json` file.
- `dataset_name`: The name of the dataset.
- `source_path`: The path to the source data file.
- `metadata`: Optional dictionary of metadata to store.
- `is_imzML`: Boolean indicating if the source is an imzML file.
"""
function update_registry(registry_path, dataset_name, source_path, metadata=nothing, is_imzML=false)
    registry = load_registry(registry_path)
    
    # Get existing entry if it exists, otherwise create new one
    existing_entry = get(registry, dataset_name, Dict{String,Any}())
    
    # Start with existing data and update only the basic fields
    entry = copy(existing_entry)
    entry["source_path"] = source_path
    entry["processed_date"] = string(now())
    entry["is_imzML"] = is_imzML
    
    # Only update metadata if provided
    if metadata !== nothing
        entry["metadata"] = metadata
    end
    
    # Note: has_mask and mask_path are preserved from existing_entry if they exist
    
    registry[dataset_name] = entry

    try
        open(registry_path, "w") do f
            JSON.print(f, registry, 4)
        end
    catch e
        @error "Failed to write to registry.json: $e"
    end
end

"""
    process_file_safely(file_path, masses, params, progress_message_ref, overall_progress_ref)

Safely processes a single MSI data file, generating and saving m/z slices.

This function handles loading data, generating slices, saving results as bitmaps,
and updating the data registry. It includes error handling and memory cleanup.

# Arguments
- `file_path`: Path to the `.imzML` file.
- `masses`: A vector of m/z values to generate slices for.
- `params`: A structure or dictionary containing processing parameters.
- `progress_message_ref`: A reference to update with progress messages.
- `overall_progress_ref`: A reference to update with overall progress.

# Returns
- A tuple `(success::Bool, message::String)`.
"""
function process_file_safely(file_path, masses, params, progress_message_ref, overall_progress_ref; use_mask::Bool=false)
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

        # --- Get mask path and load mask_matrix only if use_mask is true ---
        local mask_matrix_for_triq::Union{BitMatrix, Nothing} = nothing
        mask_path = nothing
        if use_mask
            registry = load_registry(params.registry)
            entry = get(registry, dataset_name, nothing)
            if entry !== nothing && get(entry, "has_mask", false)
                mask_path_candidate = get(entry, "mask_path", "")
                if isfile(mask_path_candidate)
                    mask_path = mask_path_candidate
                    # Load the mask matrix here
                    mask_matrix_for_triq = load_and_prepare_mask(mask_path, (local_msi_data.image_dims[1], local_msi_data.image_dims[2]))
                    println("DEBUG: Using mask: $(mask_path) with dimensions $(size(mask_matrix_for_triq))")
                else
                    @warn "Mask enabled but file not found: $(mask_path_candidate). Clearing invalid mask entry."
                    # Clear invalid mask entry
                    update_registry_mask_fields(params.registry, dataset_name, false, "")
                    mask_path = nothing
                end
            else
                @warn "Mask enabled but no valid mask entry found for: $(dataset_name)"
            end
        end

        # --- Generate Slices ---
        progress_message_ref = "Generating $(length(masses)) slices for $(dataset_name)..."
        slice_dict = get_multiple_mz_slices(local_msi_data, masses, params.tolerance, mask_path=mask_path)

        # --- Extract metadata ---
        metadata = extract_metadata(local_msi_data, file_path)

        # --- Save Slices ---
        mkpath(output_dir)

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
                sliceQuant = params.triqE ? TrIQ(slice, params.colorL, params.triqP, mask_matrix=mask_matrix_for_triq) : quantize_intensity(slice, params.colorL, mask_matrix=mask_matrix_for_triq)
                if params.medianF
                    sliceQuant = round.(UInt8, median_filter(sliceQuant))
                end
            end

            save_bitmap(joinpath(output_dir, bitmap_filename), sliceQuant, ViridisPalette)
            if !all(iszero, slice)
                generate_colorbar_image(slice, params.colorL, joinpath(output_dir, colorbar_filename); 
                    use_triq=params.triqE, triq_prob=params.triqP, mask_path=mask_path)
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

function update_registry_mask_fields(registry_path, dataset_name, has_mask, mask_path)
    registry = load_registry(registry_path)
    
    if haskey(registry, dataset_name)
        registry[dataset_name]["has_mask"] = has_mask
        registry[dataset_name]["mask_path"] = mask_path
    else
        # Create new entry if it doesn't exist
        registry[dataset_name] = Dict{String,Any}(
            "has_mask" => has_mask,
            "mask_path" => mask_path,
            "source_path" => "",
            "processed_date" => string(now()),
            "is_imzML" => false
        )
    end

    try
        open(registry_path, "w") do f
            JSON.print(f, registry, 4)
        end
    catch e
        @error "Failed to write to registry.json: $e"
    end
end

"""
    loadSurfacePlot(interfaceImg::String, mask_path::String)

Loads an image, smooths it, applies a mask, crops the data to the masked region, and creates a 3D Plotly surface plot.

# Arguments
- `interfaceImg`: Path to the image file.
- `mask_path`: Path to a mask file. Areas outside the mask will be removed and the plot axes will be cropped to the masked region.

# Returns
- `plotdata`: A vector containing the Plotly surface trace.
- `plotlayout`: The Plotly layout for the 3D plot.
"""
function loadSurfacePlot(interfaceImg::String, mask_path::String)
    # Load the image
    cleaned_img = replace(interfaceImg, r"\?.*" => "")
    cleaned_img = lstrip(cleaned_img, '/')
    var = joinpath("./public", cleaned_img)
    img = load(var)
    img_gray = Gray.(img) # Convert to grayscale
    img_array = Array(img_gray)
    elevation = Float32.(Array(img_array)) ./ 255.0 # Normalize between 0 and 1

    # Smooth the image
    sigma = 3.0
    kernel = Kernel.gaussian(sigma)
    elevation_smoothed = imfilter(elevation, kernel)

    # --- APPLY MASK and CROP DATA ---
    mask_matrix = load_and_prepare_mask(mask_path, (size(elevation_smoothed, 2), size(elevation_smoothed, 1)))
    elevation_smoothed[.!mask_matrix] .= NaN32

    non_nan_indices = findall(!isnan, elevation_smoothed)
    if isempty(non_nan_indices)
        # Return an empty plot if mask removes everything
        return [PlotlyBase.surface()], PlotlyBase.Layout(title="Empty plot: Mask covered all data")
    end

    row_indices = [idx[1] for idx in non_nan_indices]
    col_indices = [idx[2] for idx in non_nan_indices]
    
    y_range = minimum(row_indices):maximum(row_indices)
    x_range = minimum(col_indices):maximum(col_indices)

    cropped_elevation = elevation_smoothed[y_range, x_range]
    # ---

    # --- DOWNSAMPLE CROPPED DATA ---
    cropped_elevation = downsample_image(cropped_elevation, 256)
    # ---

    # Create the X, Y meshgrid coordinates for the CROPPED data
    x_coords = x_range
    y_coords = y_range
    X = repeat(reshape(x_coords, 1, length(x_coords)), length(y_coords), 1)
    Y = repeat(reshape(y_coords, length(y_coords), 1), 1, length(x_coords))

    # Define tick values and text for colorbars, ignoring NaNs
    non_nan_values = filter(!isnan, cropped_elevation)
    min_val = isempty(non_nan_values) ? 0.0 : minimum(non_nan_values)
    max_val = isempty(non_nan_values) ? 1.0 : maximum(non_nan_values)
    tickV = range(min_val, stop=max_val, length=8)
    tickT = log_tick_formatter(collect(tickV))

    # Transpose the elevation_smoothed array if Y axis is longer than X axis to fix chopping
    cropped_elevation = transpose(cropped_elevation)
    if size(cropped_elevation, 1) < size(cropped_elevation, 2)
        Y = -Y
    else
        X = -X
    end

    # Calculate the number of ticks and aspect ratio for the 3d plot
    x_nticks = min(20, length(x_coords))
    y_nticks = min(20, length(y_coords))
    z_nticks = 5
    aspect_ratio = attr(x=1, y=length(y_coords) / length(x_coords), z=0.5)
    
    # Define the layout for the 3D plot
    layout3D = PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="3D surface plot of $cleaned_img (masked, cropped, downsampled)",
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
        margin=attr(l=0, r=0, t=120, b=0, pad=0)
    )

    trace3D = PlotlyBase.surface(
        x=X[1, :],
        y=Y[:, 1],
        z=cropped_elevation,
        contours_z=attr(
            show=true,
            usecolormap=true,
            highlightcolor="limegreen",
            project_z=true
        ),
        colorscale="Viridis",
        colorbar=attr(
            tickvals=tickV,
            ticktext=tickT,
            nticks=8
        )
    )
    plotdata = [trace3D]
    plotlayout = layout3D
    return plotdata, plotlayout
end