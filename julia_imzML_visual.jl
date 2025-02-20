# IntQuantCl is originally a function of julia mzMl imzML with the adition 
# of altering the scale of the colors according to colorlevel.
function IntQuantCl( slice , colorLevel)
    # Compute scale factor for amplitude discretization
    lower  = minimum( slice )
    scale  = colorLevel / maximum( slice )
    dim    = size( slice )
    image  = zeros( UInt8, dim[1], dim[2] )
    for i in 1:length( slice )
        image[i] = convert( UInt8, floor( slice[i] * scale + 0.5 ) )
    end
    return image
end

# SaveBitmap originally a function of the mzML imzML library in julia,
# This function dinamically adjust the color palete adjusting to the ammount of colors
# available in pixmap
function SaveBitmapCl( name, pixMap::Array{UInt8,2}, colorTable::Array{UInt32,1} )
    # Get image dimensions
    dim = size( pixMap )
    if length( dim ) != 2
      return 0
    end
    # Normalize pixel values to get a more accurate reading of the image
    minVal = minimum(pixMap)
    maxVal = maximum(pixMap)
    pixMap = round.(UInt8, 255 * (pixMap .- minVal) ./ (maxVal - minVal))
    # Compute row padding
    padding = ( 4 - dim[1] & 0x3 ) & 0x3
    # Compute file dimensions. Header = 14 + 40 + ( 256 * 4 ) = 1078
    offset   = 1078
    imgBytes = dim[2] * ( dim[1] + padding )
    # Create file
    stream = open( name, "w" )
    # Save file header
    write( stream, UInt16( 0x4D42 ) )
    write( stream, UInt32[ offset + imgBytes, 0 , offset ] )
    # Save info header
    write( stream, UInt32[ 40, dim[1], dim[2], 0x80001, 0 ] )
    write( stream, UInt32[ imgBytes, 0, 0, 256, 0 ] )
    # Save color table
    write( stream, colorTable )
    if length( colorTable ) < 256
        fixTable = zeros( UInt32, 256 - length( colorTable ) )
        write( stream, fixTable )
    end
    # Save image pixels
    if padding == 0
        for i = 1:dim[2]
            write( stream, pixMap[:,i] )
        end
    else
        zeroPad = zeros( UInt8, padding )
        for i in 1:dim[2]
            write( stream, pixMap[:,i] )
            write( stream, zeroPad )
        end
    end
    # Close file
    close( stream )
end

# SaveBitmap originally a function of the mzML imzML library in julia,
# now has an adjustment for NaN values in case they exist to mantain data integrity
function GetMzSliceJl(imzML, mass, tolerance)
    # Alloc space for slice
    width = maximum(imzML[1, :])
    height = maximum(imzML[2, :])
    image = fill(0.0, width, height)

    for i in 1:size(imzML)[2]
        index = julia_mzML_imzML.FindMass(imzML[3, i], mass, tolerance)
        if index != 0
            image[imzML[1, i], imzML[2, i]] = imzML[4, i][index]
        end
    end
    # Adjustment for NaN values with 0
    replace!(image, NaN => 0.0)
    return image
end

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
            scaleanchor="y"
        ),
        yaxis=PlotlyBase.attr(
            visible=false
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
        colorbar=attr(
            tickformat=".2g"
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
        colorbar=attr(
            tickformat=".2g"
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

# Median filter: an adaptation of the R medianfilter, which averages the matrix
# with the close pixels just from the sides to reduce noise.
# This one in particular is a midpoint fiter from a 3x3 neighbour area
function medianFilterjl(pixMap)
    height, width = size(pixMap)
    target = zeros(eltype(pixMap), height, width)
    for j in 2:(width-1)
        for i in 2:(height-1)
            neighbors = []
            for dj in max(1, j-1):min(j+1, width)
                for di in max(1, i-1):min(i+1, height)
                    push!(neighbors, pixMap[di, dj])
                end
            end
            target[i, j] = median(neighbors)
        end
    end
    return target
end

# sumSpectrumPlot recieves the local directory of the image as a string,
# returns the layout and data for the surface plotly plot
# this function loads the spectra data and makes a mean to display
# its values in the spectrum plot
function sumSpectrumPlot(mzmlRoute::String)
    spectraMz=LoadMzml(mzmlRoute)
    xSpectraMz=Float64[]
    ySpectraMz=Float64[]

    layout=PlotlyBase.Layout(
        title="SUM spectrum plot",
        xaxis=PlotlyBase.attr(
            title="<i>m/z</i>",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Intensity",
            showgrid=true
        ),
        autosize=false,
        margin=attr(l=0,r=0,t=120,b=0,pad=0)
    )
    try
        xSpectraMz=mean(spectraMz[1,:])
        ySpectraMz=mean(spectraMz[2,:])
    catch e
        xSpectraMz=spectraMz[1,1]
        ySpectraMz=spectraMz[2,1]
    end
    trace=PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, mode="lines")
    plotdata=[trace] # We add the data from spectra to the plot
    plotlayout=layout
    return plotdata, plotlayout, xSpectraMz, ySpectraMz
end