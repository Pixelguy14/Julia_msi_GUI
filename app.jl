module App
# ==Packages ==
using GenieFramework # Set up Genie development environment.
using Pkg
using Libz
using PlotlyBase
using CairoMakie
using Colors
using julia_mzML_imzML
using Statistics
using NaturalSort
using Images
using LinearAlgebra
using NativeFileDialog # Opens the file explorer depending on the OS
using StipplePlotly
using Printf
include("./julia_imzML_visual.jl")
@genietools

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded

# == Search functions ==
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
    #println("type of img: $(typeof(img))")
    # Convert to grayscale
    img_gray=Gray.(img)
    img_array=Array(img_gray)
    #println(typeof(img_array))
    elevation=Float32.(Array(img_array)) ./ 255.0
    #println("type of elevation: $(typeof(elevation))")
    # Get the X, Y coordinates of the image
    height, width=size(img_array)
    #println("height: $(height), width: $(width)")
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
        )
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
    # parse matrix to int32 to closely resemble similarity to original bitmap
    elevation=int32.(Array(img_gray))  
    
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
        title="2D Topographic Map",
        xaxis=PlotlyBase.attr(
            title="X",
            scaleanchor="y"
        ),
        yaxis=PlotlyBase.attr(
            title="Y"
        ),
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
    #println("Image type:", typeof(img))
    img_gray=Gray.(img) # Convert to grayscale
    #println("Grayscale image type:", typeof(img_gray)) 
    img_array=Array(img_gray)
    elevation=Float32.(Array(img_gray)) ./ 255.0 # Normalize between 0 and 1
    #println("Elevation size:", size(elevation)) 
    # Smooth the image 
    sigma=3.0
    kernel=Kernel.gaussian(sigma)
    #println(size(kernel))
    elevation_smoothed=imfilter(elevation, kernel) 
    #println("Smoothed elevation size:", size(elevation_smoothed))
    # Transpose the elevation_smoothed array
    # Create the X, Y meshgrid coordinates 
    x=1:size(elevation_smoothed, 2) 
    y=1:size(elevation_smoothed, 1)
    X=repeat(reshape(x, 1, length(x)), length(y), 1)
    #println("Size of X:", size(X))
    Y=repeat(reshape(y, length(y), 1), 1, length(x))
    #println("Size of Y:", size(Y))
    # Calculate the number of ticks and aspect ratio for the 3d plot
    x_nticks=min(20, length(x))
    y_nticks=min(20, length(y))
    z_nticks=5
    aspect_ratio=attr(x=1, y=length(y) / length(x), z=0.5)
    # Define the layout for the 3D plot
    layout3D=PlotlyBase.Layout(
        title="3D Surface Plot",
        scene=attr(
            xaxis_nticks=x_nticks,
            yaxis_nticks=y_nticks,
            zaxis_nticks=z_nticks, 
            camera=attr(eye=attr(x=0, y=-1, z=0.5)), 
            aspectratio=aspect_ratio
        )
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

# == Reactive code ==
# Reactive code to make the UI interactive
@app begin
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block
    #@out test="/test.bmp" #slash means it's getting the info from 'public' folder

    ## Interface non Variables
    @out btnStartDisable=true
    @out btnPlotDisable=false
    @out btnSpectraDisable=true
    # Loading animations
    @in progress=false
    @in progressPlot=false
    @in progressSpectraPlot=false
    # Text field validations
    @in triqEnabled=false
    @in SpectraEnabled=false
    # Dialogs
    @in warning_msg=false
    @in CompareDialog=false

    ## Interface Variables
    @in file_route=""
    @in file_name=""
    @in Nmass=0.0
    @in Tol=0.1
    @in triqProb=0.98
    @in colorLevel=20

    ## Interface Buttons
    @in btnSearch=false # To search for files in your device
    @in mainProcess=false # To generate images
    @in compareBtn=false # To open dialog
    @in createSumPlot=false # To generate sum spectrum plot
    @in createXYPlot=false # To generate an spectrum plot according to the xy values inputed
    @in image3dPlot=false # To generate 3d plot based on current image
    @in triq3dPlot=false # To generate 3d plot based on current triq image
    @in imageCPlot=false # To generate contour plots of current image
    @in triqCPlot=false # To generate contour plots of current triq image
    # Image change buttons
    @in imgPlus=false
    @in imgMinus=false
    @in imgPlusT=false
    @in imgMinusT=false

    ## Tabulation variables
    @out tabIDs=["tab0","tab1","tab2","tab3","tab4"]
    @out tabLabels=["Image", "TrIQ", "Spectrum Plot", "Topology Plot","Surface Plot"]
    @in selectedTab="tab0"
    @out CompTabIDs=["tab0","tab1","tab2","tab3","tab4"]
    @out CompTabLabels=["Image", "TrIQ", "Spectrum Plot", "Topology Plot","Surface Plot"]
    @in CompSelectedTab="tab0"

    # Interface Images
    @out imgInt="/.bmp" # image Interface
    @out imgIntT="/.bmp" # image Interface TrIQ
    @out colorbar="/.png"
    @out colorbarT="/.png"
    @out imgWidth=0
    @out imgHeight=0

    # Messages to interface variables
    @out msg=""
    @out msgimg=""
    @out msgtriq=""

    # Saves the route where imzML and mzML files are located
    @out full_route=""
    @out full_routeMz=""
    @out full_routeMz2=""

    # For the creation of images with a more specific mass charge
    @out text_nmass=""

    # For image search image lists we apply a filter that searches specific type of images into our public folder, then we sort it in a "numerical" order
    @in msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    @in triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    
    # Set current image for the list to display
    @out current_msi=""
    @out current_col_msi=""
    @out current_triq=""
    @out current_col_triq=""

    ## Time measurement variables
    @out sTime=time()
    @out fTime=time()
    @out eTime=time()

    ## Plots
    # Local image to plot
    layoutImg=PlotlyBase.Layout(
        xaxis=PlotlyBase.attr(
            visible=false,
            scaleanchor="y"
        ),
        yaxis=PlotlyBase.attr(
            visible=false
        )
    )
    traceImg=PlotlyBase.heatmap(x=[], y=[])
    @out plotdataImg=[traceImg]
    @out plotlayoutImg=layoutImg
    # For triq image 
    @out plotdataImgT=[traceImg]
    @out plotlayoutImgT=layoutImg
    # Interface Plot Spectrum
    layoutSpectra=PlotlyBase.Layout(
        title="SUM Spectrum plot",
        xaxis=PlotlyBase.attr(
            title="<i>m/z</i>",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Intensity",
            showgrid=true
        )
    )
    # Dummy 2D scatter plot
    traceSpectra=PlotlyBase.scatter(x=[], y=[], mode="lines")
    # Create conection to frontend
    @out plotdata=[traceSpectra]
    @out plotlayout=layoutSpectra
    @in xCoord=0
    @in yCoord=0
    @out xSpectraMz=Float64[]
    @out ySpectraMz=Float64[]

    # Interactive plot reactions
    @in data_click=Dict{String,Any}()
    #@in data_selected=Dict{String,Any}() # Selected is for areas, this can work for the masks
    #<plotly id="plotStyle" :data="plotdata" :layout="plotlayout" @click="data_selected" class="q-pa-none q-ma-none sync_data"></plotly>

    # Interface Plot Surface
    layoutContour=PlotlyBase.Layout(
        title="2D Topographic Map",
        xaxis=PlotlyBase.attr(
            title="X",
            scaleanchor="y"
            ),
            yaxis=PlotlyBase.attr(
                title="Y"
                ),
        )
    # Dummy 2D surface plot
    traceContour=PlotlyBase.scatter(x=[], y=[], mode="lines")
    # Create conection to frontend
    @out plotdataC=[traceContour]
    @out plotlayoutC=layoutContour

    # Interface Plot 3d
    # Define the layout for the 3D plot
    layout3D=PlotlyBase.Layout(
        title="3D Surface Plot",
        scene=attr(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            xaxis_nticks=20, 
            yaxis_nticks=20, 
            zaxis_nticks=4, 
            camera=attr(eye=attr(x=0, y=-1, z=0.5)), 
            aspectratio=attr(x=1, y=1, z=0.2)
        )
    )

    # Dummy 3D surface plot
    x=1:10
    y=1:10
    z=[sin(i * j / 10) for i in x, j in y]
    trace3D=PlotlyBase.surface(x=[], y=[], z=[],
                                contours_z=attr(
                                    show=true,
                                    usecolormap=true,
                                    highlightcolor="limegreen",
                                    project_z=true
                                ), colorscale="Viridis")
    # Create conection to frontend
    @out plotdata3d=[trace3D]
    @out plotlayout3d=layout3D

    # == Reactive handlers ==
    # Reactive handlers watch a variable and execute a block of code when its value changes
    # The onbutton handler will set the variable to false after the block is executed

    @onbutton btnSearch begin
        full_route=pick_file(; filterlist="imzML,mzML")
        if isnothing(full_route)
            #println("No file selected")
            msg="No file selected"
            warning_msg=true
            btnStartDisable=true
        else
            if endswith(full_route, "imzML") # Case if the file loaded is imzML
                #println("Selected file path: ", full_route)
                btnStartDisable=false
                btnPlotDisable=false
                # Splitting the route with regex from imzml to mzml so the plotting can work
                full_routeMz=replace(full_route, r"\.[^.]*$" => ".mzML") 
                if isfile(full_routeMz)
                    # We enable coord search and spectra plot creation
                    btnSpectraDisable=false
                    SpectraEnabled=true
                else
                    # If there's no MzML file, we deny access again
                    btnSpectraDisable=true
                    SpectraEnabled=false
                end
            else # Case if the file loaded is mzML
                full_routeMz=full_route
                btnSpectraDisable=false
                SpectraEnabled=true
                # Splitting the route the same way
                full_route=replace(full_route, r"\.[^.]*$" => ".imzML") 
                if isfile(full_route)
                    btnStartDisable=false
                else
                    btnStartDisable=true
                    full_route=full_routeMz
                end
            end
            xCoord=0
            yCoord=0
        end
    end
    
    @onbutton mainProcess begin
    #@onchange Nmass begin
        progress=true # Start progress button animation
        btnStartDisable=true # We disable the button to avoid multiple requests
        btnPlotDisable=true
        btnSpectraDisable=true
        text_nmass=replace(string(Nmass), "." => "_")
        sTime=time()
        #full_route=joinpath(file_route, file_name)
        if isfile(full_route) && Nmass > 0 && Tol > 0 && Tol <=1 && colorLevel > 1 && colorLevel < 257
            msg="File exists, Nmass=$(Nmass) Tol=$(Tol). Loading file will begin, please be patient."
            try
                spectra=LoadImzml(full_route)
                msg="File loaded. Creating Spectra with the specific mass and tolerance, please be patient."
                slice=GetSlice(spectra, Nmass, Tol)
                fig=CairoMakie.Figure(size=(128, 256)) # Container
                # Append a query string to force the image to refresh 
                timestamp=string(time_ns()) 
                if triqEnabled # If we have TrIQ
                    if triqProb < 0.8 || triqProb > 1
                        msg="Incorrect TrIQ values, please adjust accordingly and try again."
                        warning_msg=true
                    else
                        image_path=joinpath("./public", "TrIQ_$(text_nmass).bmp")
                        sliceTriq=TrIQ(slice, colorLevel, triqProb)
                        #println("slice raw: $(typeof(slice))")
                        #println("TriQ matrix: $(typeof(sliceTriq))")
                        sliceTriq=reverse(sliceTriq, dims=2)
                        ##SaveBitmap(joinpath("public", "TrIQ_$(text_nmass).bmp"),sliceTriq,ViridisPalette)
                        SaveBitmapCl(joinpath("public", "TrIQ_$(text_nmass).bmp"),sliceTriq,ViridisPalette)
                        # Use timestamp to refresh image interface container
                        imgIntT="/TrIQ_$(text_nmass).bmp?t=$(timestamp)"
                        plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
                        # Get current image 
                        current_triq="TrIQ_$(text_nmass).bmp"
                        msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
                        # Create colorbar
                        #ticks=round.(range(0, (stop=maximum(slice)*triqProb), length=10), sigdigits=3)
                        ticks=range(0, (stop=maximum(slice)*triqProb), length=5)

                        """
                        function log_tick_formatter(values)
                            return map(v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)), values)
                        end
                        """

                        println("ticks: $(typeof(ticks))")
                        #Colorbar(fig[1, 1], colormap=cgrad(:viridis, colorLevel, categorical=true), limits=(0, maximum(slice)*triqProb),ticks=ticks,tickformat=log_tick_formatter, label="Intensity", size = 25)
                        Colorbar(fig[1, 1], colormap=cgrad(:viridis, colorLevel, categorical=true), limits=(0, maximum(slice)*triqProb),ticks=ticks, label="Intensity", size = 25)
                        println("Colorbar: $(typeof(fig))")
                        save("public/colorbar_TrIQ_$(text_nmass).png", fig)
                        colorbarT="/colorbar_TrIQ_$(text_nmass).png?t=$(timestamp)"
                        # Get current colorbar 
                        current_col_triq="colorbar_TrIQ_$(text_nmass).png"
                        # We update the directory to include the new placed images.
                        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
                        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
                        fTime=time()
                        eTime=round(fTime-sTime,digits=3)
                        msg="The file has been created in $(eTime) seconds successfully inside the 'public' folder of the app"
                        selectedTab="tab1"
                        #println("all msi in folder=",triq_bmp)
                        #println("all col msi in folder=",col_triq_png)
                    end
                else # If we don't use TrIQ
                    image_path=joinpath("./public", "MSI_$(text_nmass).bmp")
                    ##sliceQuant=IntQuant(slice)
                    sliceQuant=IntQuantCl(slice,Int(colorLevel-1))
                    #println("slice raw: $(typeof(slice))")
                    #println("slice in intQuant: $(typeof(sliceQuant))")
                    sliceQuant=reverse(sliceQuant, dims=2)
                    ##SaveBitmap(joinpath("public", "MSI_$(text_nmass).bmp"),sliceQuant,ViridisPalette)
                    SaveBitmapCl(joinpath("public", "MSI_$(text_nmass).bmp"),sliceQuant,ViridisPalette)
                    # Use timestamp to refresh image interface container
                    imgInt="/MSI_$(text_nmass).bmp?t=$(timestamp)"
                    plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
                    # Get current image 
                    current_msi="MSI_$(text_nmass).bmp"
                    msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
                    # Create colorbar 
                    #ticks=round.(range(0, stop=maximum(slice), length=10), sigdigits=3)
                    ticks=range(0, stop=maximum(slice), length=5)
                    Colorbar(fig[1, 1], colormap=cgrad(:viridis, colorLevel, categorical=true), limits=(0, maximum(slice)),ticks=ticks, label="Intensity", size = 25)
                    save("public/colorbar_MSI_$(text_nmass).png", fig)
                    colorbar="/colorbar_MSI_$(text_nmass).png?t=$(timestamp)"
                    # Get current colorbar 
                    current_col_msi="colorbar_MSI_$(text_nmass).png"
                    # We update the directory to include the new placed images.
                    msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
                    col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
                    selectedTab="tab0"
                    fTime=time()
                    eTime=round(fTime-sTime,digits=3)
                    msg="The file has been created in $(eTime) seconds successfully inside the 'public' folder of the app"
                    #println("all msi in folder=",msi_bmp)
                    #println("all col msi in folder=",col_msi_png)
                end
            catch e
                msg="There was an error loading the ImzML file, please verify the file accordingly and try again. $(e)"
                warning_msg=true
            end
        else
            msg="File does not exist or a parameter is incorrect, please try again."
            warning_msg=true
        end
        spectra=nothing # Important for memory cleaning
        slice=nothing
        GC.gc() # Trigger garbage collection
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
        end
        btnStartDisable=false
        btnPlotDisable=false
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
        progress=false
    end

    @onbutton createSumPlot begin
        msg="Sum spectrum plot selected"
        sTime=time()
        #full_route=joinpath( file_route, file_name )
        if isfile(full_routeMz) # Check if the file exists
            progressSpectraPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            msg="Loading plot..."
            spectraMz=LoadMzml(full_routeMz)
            layoutSpectra=PlotlyBase.Layout(
                title="SUM Spectrum plot",
                xaxis=PlotlyBase.attr(
                    title="<i>m/z</i>",
                    showgrid=true
                ),
                yaxis=PlotlyBase.attr(
                    title="Intensity",
                    showgrid=true
                ),
                autosize=false
            )
            # dims=size(spectraMz)
            # scansMax=dims[2] # we get the total of scansMax
            xSpectraMz=mean(spectraMz[1,:])
            ySpectraMz=mean(spectraMz[2,:])
            traceSpectra=PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, mode="lines")
            plotdata=[traceSpectra] # We add the data from spectra to the plot
            plotlayout=layoutSpectra
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end
            selectedTab="tab2"
            fTime=time()
            eTime=round(fTime-sTime,digits=3)
            msg="Plot loaded in $(eTime) seconds"
        else
            msg="there was an error with the mzML, please try again"
            warning_msg=true
        end
        progressSpectraPlot=false
        btnPlotDisable=false
        if endswith(full_route, "imzML")
            btnStartDisable=false
        end
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
    end

    @onbutton createXYPlot begin
        msg="Sum spectrum plot selected"
        sTime=time()
        #full_route=joinpath( file_route, file_name )
        if isfile(full_routeMz) # Check if the file exists
            progressSpectraPlot=true
            btnStartDisable=true
            btnPlotDisable=true
            btnSpectraDisable=true
            msg="Loading plot..."
            spectraMz=LoadMzml(full_routeMz)
            layoutSpectra=PlotlyBase.Layout(
                title="($xCoord, $yCoord) Specific spectrum plot",
                xaxis=PlotlyBase.attr(
                    title="<i>m/z</i>",
                    showgrid=true
                ),
                yaxis=PlotlyBase.attr(
                    title="Intensity",
                    showgrid=true
                ),
                autosize=false
            )
            xSpectraMz=spectraMz[1,abs(xCoord)]
            ySpectraMz=spectraMz[2,abs(yCoord)]
            traceSpectra=PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, mode="lines")
            plotdata=[traceSpectra] # We add the data from spectra to the plot
            plotlayout=layoutSpectra
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end
            selectedTab="tab2"
            fTime=time()
            eTime=round(fTime-sTime,digits=3)
            msg="Plot loaded in $(eTime) seconds"
        else
            msg="there was an error with the mzML or the coordenates, please try again"
            warning_msg=true
        end
        progressSpectraPlot=false
        btnPlotDisable=false
        if endswith(full_route, "imzML")
            btnStartDisable=false
        end
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
    end

    # Image loaders based on the position of the current image (increment and decrement for both normal and filter)
    # And a pre-generated list from all image files from /public folder
    @onbutton imgMinus begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images listed in the public folder
        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        
        new_msi=decrement_image(current_msi, msi_bmp)
        new_col_msi=decrement_image(current_col_msi, col_msi_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_msi=new_msi
            current_col_msi=new_col_msi
            imgInt="/$(current_msi)?t=$(timestamp)"
            colorbar="/$(current_col_msi)?t=$(timestamp)"

            text_nmass=replace(current_msi, "MSI_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImg=[traceImg]
            msgimg=""
        end
    end
    @onbutton imgPlus begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns())
        # Update the array of images listed in the public folder
        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        
        new_msi=increment_image(current_msi, msi_bmp)
        new_col_msi=increment_image(current_col_msi, col_msi_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_msi=new_msi
            current_col_msi=new_col_msi
            imgInt="/$(current_msi)?t=$(timestamp)"
            colorbar="/$(current_col_msi)?t=$(timestamp)"

            text_nmass=replace(current_msi, "MSI_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImg=[traceImg]
            msgimg=""
        end
    end

    @onbutton imgMinusT begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images with TrIQ filter listed in the public folder
        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)

        new_msi=decrement_image(current_triq, triq_bmp)
        new_col_msi=decrement_image(current_col_triq, col_triq_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_triq=new_msi
            current_col_triq=new_col_msi
            imgIntT="/$(current_triq)?t=$(timestamp)"
            colorbarT="/$(current_col_triq)?t=$(timestamp)"

            text_nmass=replace(current_triq, "TrIQ_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgT=[traceImg]
            msgtriq=""
        end
    end
    @onbutton imgPlusT begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images with TrIQ filter listed in the public folder
        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)

        new_msi=increment_image(current_triq, triq_bmp)
        new_col_msi=increment_image(current_col_triq, col_triq_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_triq=new_msi
            current_col_triq=new_col_msi
            imgIntT="/$(current_triq)?t=$(timestamp)"
            colorbarT="/$(current_col_triq)?t=$(timestamp)"

            text_nmass=replace(current_triq, "TrIQ_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgT=[traceImg]
            msgtriq=""
        end
    end

    # 3d plot
    @onbutton image3dPlot begin
        msg="Image 3D plot selected"
        cleaned_imgInt=replace(imgInt, r"\?.*" => "")
        cleaned_imgInt=lstrip(cleaned_imgInt, '/')
        var=joinpath( "./public", cleaned_imgInt )
        sTime=time()

        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            btnSpectraDisable=true
            try
                plotdata3d, plotlayout3d=loadSurfacePlot(imgInt)
                GC.gc() # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                end
                selectedTab="tab4"
                fTime=time()
                eTime=round(fTime-sTime,digits=3)
                msg="Plot loaded in $(eTime) seconds"
            catch e 
                msg="Failed to load and process image: $e" 
                warning_msg=true 
                println(msg)
            end
        else
            msg="Image could not be 3d plotted"
            warning_msg=true
        end
        progressPlot=false
        btnPlotDisable=false
        btnStartDisable=false
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
    end
    # 3d plot for TrIQ
    @onbutton triq3dPlot begin
        msg="TrIQ 3D plot selected"
        cleaned_imgIntT=replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT=lstrip(cleaned_imgIntT, '/')
        var=joinpath( "./public", cleaned_imgIntT )
        sTime=time()

        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            btnSpectraDisable=true
            try
                plotdata3d, plotlayout3d=loadSurfacePlot(imgIntT)
                GC.gc() # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                end
                selectedTab="tab4"
                fTime=time()
                eTime=round(fTime-sTime,digits=3)
                msg="Plot loaded in $(eTime) seconds"
            catch e 
                msg="Failed to load and process image: $e" 
                warning_msg=true 
                println(msg)
            end
        else
            msg="Image could not be 3d plotted"
            warning_msg=true
        end
        progressPlot=false
        btnPlotDisable=false
        btnStartDisable=false
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
    end

    # Contour 2d plot
    @onbutton imageCPlot begin
        msg="Image 2D plot selected"
        cleaned_imgInt=replace(imgInt, r"\?.*" => "")
        cleaned_imgInt=lstrip(cleaned_imgInt, '/')
        var=joinpath("./public", cleaned_imgInt)
        sTime=time()
    
        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            btnSpectraDisable=true
            try
                img=load(var)
                plotdataC,plotlayoutC=loadContourPlot(imgInt)
                GC.gc()  # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0)  # Ensure Julia returns the freed memory to OS
                end
                selectedTab="tab3"
                fTime=time()
                eTime=round(fTime-sTime,digits=3)
                msg="Plot loaded in $(eTime) seconds"
            catch e
                msg="Failed to load and process image: $e"
                warning_msg=true
                println(msg)
            end
        else
            msg="Image could not be 2D plotted"
            warning_msg=true
        end
        progressPlot=false
        btnPlotDisable=false
        btnStartDisable=false
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
    end
    # Contour 2d plot for TrIQ
    @onbutton triqCPlot begin
        msg="Image 2D plot selected"
        cleaned_imgIntT=replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT=lstrip(cleaned_imgIntT, '/')
        var=joinpath("./public", cleaned_imgIntT)
        sTime=time()
    
        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            btnSpectraDisable=true
            try
                img=load(var)
                plotdataC,plotlayoutC=loadContourPlot(imgIntT)
                GC.gc()  # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0)  # Ensure Julia returns the freed memory to OS
                end
                selectedTab="tab3"
                fTime=time()
                eTime=round(fTime-sTime,digits=3)
                msg="Plot loaded in $(eTime) seconds"
            catch e
                msg="Failed to load and process image: $e"
                warning_msg=true
                println(msg)
            end
        else
            msg="Image could not be 2D plotted"
            warning_msg=true
        end
        progressPlot=false
        btnPlotDisable=false
        btnStartDisable=false
        if isfile(full_routeMz)
            # We enable coord search and spectra plot creation
            btnSpectraDisable=false
            SpectraEnabled=true
        end
    end

    @onbutton compareBtn begin
        CompareDialog=true
        # We remove the red lines
        traceSpectra=PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, mode="lines",name="Spectra",showlegend=false)
        plotdata=[traceSpectra]
    end

    # Event detection for clicking on the spectrum plot
    @onchange data_click begin
        if selectedTab == "tab2"
            if !isempty(xSpectraMz)
                #println("Clicked data on sum spectrum plot : ", data_click)
                spectracoords=reshape(plotdata, 1, length(plotdata))
                #println("Spectra: $(ndims(spectracoords))")
                # Extract x and y values from data_click 
                cursor_data=data_click["cursor"] 
                x_value=cursor_data["x"]
                y_value=cursor_data["y"]  # Get the x and y values from the click of the cursor
                closest_distance=Inf
                
                for val in spectracoords
                    # Find the index where x is within a range
                    start_idx=findfirst(x -> x >= x_value - 10, val[:x])
                    end_idx=findlast(x -> x <= x_value + 10, val[:x])
                    
                    # Ensure the index are valid and within range
                    if start_idx !== nothing && end_idx !== nothing
                        for i in start_idx:end_idx
                            spectra_x=val[:x][i]
                            spectra_y=val[:y][i]
                            distance=sqrt((spectra_x - x_value)^2 + (spectra_y - y_value)^2)  # Calculate distance
                            if distance < closest_distance
                                closest_distance=distance
                                Nmass=round(spectra_x, digits=2)
                            end
                        end
                    end
                end
                layoutSpectra=PlotlyBase.Layout(
                            title="SUM Spectrum plot",
                            xaxis=PlotlyBase.attr(
                                title="<i>m/z</i>",
                                showgrid=true
                            ),
                            yaxis=PlotlyBase.attr(
                                title="Intensity",
                                showgrid=true
                            ),
                            autosize=false
                    )
                traceSpectra=PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, mode="lines",name="Spectra",showlegend=false)
                trace2=PlotlyBase.scatter(x=[Nmass, Nmass],y=[0, maximum(ySpectraMz)],mode="lines",line=attr(color="red", width=0.5),name="<i>m/z</i> selected",showlegend=false)
                plotdata=[traceSpectra,trace2] # We add the data from spectra and the red line to the plot
                plotlayout=layoutSpectra
            end
        elseif selectedTab == "tab1"
            #println("you have clicked the triq image")
            cursor_data=data_click["cursor"] 
            xCoord=Int32(round(cursor_data["x"]))
            if xCoord < 0
                xCoord=0
            elseif xCoord > imgWidth
                xCoord=imgWidth
            end
            yCoord=Int32(round(cursor_data["y"]))
            if yCoord > 0
                yCoord=0
            elseif yCoord < -imgHeight
                yCoord=-imgHeight
            end # Get the x and y values from the click of the cursor and make sure they don't exceed image proportions
            #plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
            plotdataImgT=filter(trace -> !(get(trace, :name, "") in ["Line X", "Line Y"]), plotdataImgT)
            trace1, trace2=crossLinesPlot(xCoord, yCoord, imgWidth, -imgHeight)
            plotdataImgT=append!(plotdataImgT, [trace1, trace2])
        elseif selectedTab == "tab0"
            #println("you have clicked the normal image")
            cursor_data=data_click["cursor"] 
            xCoord=Int32(round(cursor_data["x"]))
            if xCoord < 0
                xCoord=0
            elseif xCoord > imgWidth
                xCoord=imgWidth
            end
            yCoord=Int32(round(cursor_data["y"]))  
            if yCoord > 0
                yCoord=0
            elseif yCoord < -imgHeight
                yCoord=-imgHeight
            end # Get the x and y values from the click of the cursor and make sure they don't exceed image proportions
            #plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
            plotdataImg=filter(trace -> !(get(trace, :name, "") in ["Line X", "Line Y"]), plotdataImg)
            trace1, trace2=crossLinesPlot(xCoord, yCoord, imgWidth, -imgHeight)
            plotdataImg=append!(plotdataImg, [trace1, trace2])
        end
    end

    # WIP add an x and y input and a plot to select pixels in an image and then calculate a plot (not the sum of plots)
    # optional: red lines to signal which pixel is selected in the image plot

    @mounted watchplots()

    GC.gc() # Trigger garbage collection
    if Sys.islinux()
        ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
    end
end
# == Pages ==
# Register a new route and the page that will be loaded on access
@page("/", "app.jl.html")
end

# == Advanced features ==
#=
- The @private macro defines a reactive variable that is not sent to the browser.
This is useful for storing data that is unique to each user session but is not needed
in the UI.
    @private table=DataFrame(a=1:10, b=10:19, c=20:29)

=#

