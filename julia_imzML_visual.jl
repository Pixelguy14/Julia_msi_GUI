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
        title="2D topographic map of $cleaned_img",
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
        title="3D surface plot of $cleaned_img",
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

# meanSpectrumPlot recieves the local directory of the image as a string,
# returns the layout and data for the surface plotly plot
# this function loads the spectra data and makes a mean to display
# its values in the spectrum plot
function meanSpectrumPlot(data::MSIData)
    layout = PlotlyBase.Layout(
        title="Average Spectrum Plot",
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
        trace = PlotlyBase.stem(x=Float64[], y=Float64[])
    else
        trace = PlotlyBase.stem(x=xSpectraMz, y=ySpectraMz, marker=attr(size=1, color="blue", opacity=0.5), name="Average", hoverinfo="x",hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    end


    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end

function xySpectrumPlot(data::MSIData, xCoord::Int, yCoord::Int, imgWidth::Int, imgHeight::Int)
    local mz::AbstractVector, intensity::AbstractVector
    local plot_title::String

    is_imaging = data.source isa ImzMLSource

    if is_imaging
        # For imaging data, use (X, Y) coordinates
        x = clamp(xCoord, 1, imgWidth)
        y = clamp(yCoord, 1, imgHeight)
        
        mz, intensity = GetSpectrum(data, x, y)
        plot_title = "Spectrum at ($x, $y)"
    else
        # For non-imaging data, treat xCoord as the spectrum index
        index = clamp(xCoord, 1, length(data.spectra_metadata))
        
        mz, intensity = GetSpectrum(data, index)
        plot_title = "Spectrum #$index"
    end

    layout = PlotlyBase.Layout(
        title=plot_title,
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
    
    # Return the full data for other uses, and the plot data
    return plotdata, plotlayout, mz, intensity
end

function sumSpectrumPlot(data::MSIData)
    layout = PlotlyBase.Layout(
        title="Total Spectrum Plot",
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
        trace = PlotlyBase.stem(x=Float64[], y=Float64[])
    else
        trace = PlotlyBase.stem(x=xSpectraMz, y=ySpectraMz, marker=attr(size=1, color="blue", opacity=0.5), name="Total", hoverinfo="x",hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>")
    end

    plotdata = [trace]
    plotlayout = layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end
