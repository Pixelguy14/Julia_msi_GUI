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
include("./julia_imzML_visual.jl")
@genietools

# == Reactive code ==
# Reactive code to make the UI interactive
@app begin
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block

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
    @in MFilterEnabled=false
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
    # Image change comparative buttons
    @in imgPlusComp=false
    @in imgMinusComp=false
    @in imgPlusTComp=false
    @in imgMinusTComp=false

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
    # interface controlling for the comparative view
    @out imgIntComp="/.bmp" # image Interface
    @out imgIntTComp="/.bmp" # image Interface TrIQ
    @out colorbarComp="/.png"
    @out colorbarTComp="/.png"
    @out imgWidth=0
    @out imgHeight=0

    # Optical Image Overlay & Transparency
    @in imgTrans=1.0
    @in progressOptical=false
    @out btnOpticalDisable=true
    @in btnOptical=false
    @in btnOpticalT=false
    @in opticalOverTriq=false
    @out imgRoute=""

    # Messages to interface variables
    @out msg=""
    @out msgimg=""
    @out msgtriq=""
    # Reiteration of the messages under the image to know which spectra is being visualized
    @out msgimgComp=""
    @out msgtriqComp=""

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
    # We reiterate the process to display in the comparative view
    @out current_msiComp=""
    @out current_col_msiComp=""
    @out current_triqComp=""
    @out current_col_triqComp=""

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
        ),
        margin=attr(l=0,r=0,t=0,b=0,pad=0)
    )
    traceImg=PlotlyBase.heatmap(x=[], y=[])
    @out plotdataImg=[traceImg]
    @out plotlayoutImg=layoutImg
    # For the image in the comparative view
    @out plotdataImgComp=[traceImg]
    @out plotlayoutImgComp=layoutImg

    # For triq image 
    @out plotdataImgT=[traceImg]
    @out plotlayoutImgT=layoutImg
    # For the triq image in the comparative view
    @out plotdataImgTComp=[traceImg]
    @out plotlayoutImgTComp=layoutImg
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
        ),
        margin=attr(l=0,r=0,t=120,b=0,pad=0)
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
        margin=attr(l=0,r=0,t=120,b=0,pad=0)
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
        ),
        margin=attr(l=0,r=0,t=120,b=0,pad=0)
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
            msg="No file selected"
            warning_msg=true
            btnStartDisable=true
        else
            if endswith(full_route, "imzML") # Case if the file loaded is imzML
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
        progress=true # Start progress button animation
        btnStartDisable=true # We disable the button to avoid multiple requests
        btnPlotDisable=true
        btnSpectraDisable=true
        text_nmass=replace(string(Nmass), "." => "_")
        sTime=time()
        if isfile(full_route) && Nmass > 0 && Tol > 0 && Tol <=1 && colorLevel > 1 && colorLevel < 257
            msg="File exists, Nmass=$(Nmass) Tol=$(Tol). Loading file will begin, please be patient."
            try
                spectra=LoadImzml(full_route)
                msg="File loaded. Creating spectra with the specific mass and tolerance, please be patient."
                slice=GetMzSliceJl(spectra,Nmass,Tol)
                fig=CairoMakie.Figure(size=(150, 250)) # Container
                # Append a query string to force the image to refresh 
                timestamp=string(time_ns()) 
                if triqEnabled # If we have TrIQ
                    if triqProb < 0.8 || triqProb > 1
                        msg="Incorrect TrIQ values, please adjust accordingly and try again."
                        warning_msg=true
                    else
                        image_path=joinpath("./public", "TrIQ_$(text_nmass).bmp")
                        valid_slice=false
                        while Tol <= 1.0 && !valid_slice
                            try
                                slice=GetMzSliceJl(spectra, Nmass, Tol)
                                sliceTriq=TrIQ(slice, colorLevel, triqProb)
                                if MFilterEnabled # If the Median filter is ON
                                    sliceTriq=medianFilterjl(sliceTriq)
                                end
                                valid_slice=true
                            catch e
                                msg="Warning: insufficient tolerance, inputs modified to allow the creation of an image regardless=$Tol: $e"
                                Tol += 0.1
                            end
                        end
                        sliceTriq=reverse(sliceTriq, dims=2)
                        SaveBitmapCl(joinpath("public", "TrIQ_$(text_nmass).bmp"),sliceTriq,ViridisPalette)
                        # Use timestamp to refresh image interface container
                        imgIntT="/TrIQ_$(text_nmass).bmp?t=$(timestamp)"
                        plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
                        # Get current image 
                        current_triq="TrIQ_$(text_nmass).bmp"
                        msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
                        # Create colorbar
                        bound =julia_mzML_imzML.GetOutlierThres(slice, triqProb)
                        levels=range(bound[1],stop=bound[2], length=8)
                        levels=vcat(levels, 2*levels[end]-levels[end-1])
                        Colorbar(fig[1, 1], colormap=cgrad(:viridis, colorLevel, categorical=true), limits=(0, bound[2]),ticks=levels,tickformat=log_tick_formatter, label="Intensity", size=25)
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
                    end
                else # If we don't use TrIQ
                    image_path=joinpath("./public", "MSI_$(text_nmass).bmp")
                    try
                        sliceQuant=IntQuantCl(slice,Int(colorLevel-1))
                        if MFilterEnabled # If the Median filter is ON
                            sliceQuant=medianFilterjl(sliceQuant)
                        end
                    catch e
                        msg="Warning: $e"
                    end
                    sliceQuant=reverse(sliceQuant, dims=2)
                    SaveBitmapCl(joinpath("public", "MSI_$(text_nmass).bmp"),sliceQuant,ViridisPalette)
                    # Use timestamp to refresh image interface container
                    imgInt="/MSI_$(text_nmass).bmp?t=$(timestamp)"
                    plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
                    # Get current image 
                    current_msi="MSI_$(text_nmass).bmp"
                    msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
                    # Create colorbar 
                    levels=range(0,maximum(slice),length=8)
                    Colorbar(fig[1, 1], colormap=cgrad(:viridis, colorLevel, categorical=true), limits=(0, maximum(slice)),ticks=levels,tickformat=log_tick_formatter, label="Intensity", size=25)
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
        btnOpticalDisable=false
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
                autosize=false,
                margin=attr(l=0,r=0,t=120,b=0,pad=0)
            )
            if xCoord < 1
                xCoord=1
            elseif xCoord > imgWidth
                xCoord=imgWidth
            end
            if yCoord > -1
                yCoord=-1
            elseif yCoord < -imgHeight
                yCoord=-imgHeight
            end
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
            btnOpticalDisable=false
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
            btnOpticalDisable=false
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
            btnOpticalDisable=false
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
            btnOpticalDisable=false
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgT=[traceImg]
            msgtriq=""
        end
    end

    # Image loaders for the comparative view
    @onbutton imgMinusComp begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images listed in the public folder
        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    
        new_msi=decrement_image(current_msiComp, msi_bmp)
        new_col_msi=decrement_image(current_col_msiComp, col_msi_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_msiComp=new_msi
            current_col_msiComp=new_col_msi
            imgIntComp="/$(current_msiComp)?t=$(timestamp)"
            colorbarComp="/$(current_col_msiComp)?t=$(timestamp)"
    
            text_nmass=replace(current_msiComp, "MSI_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgimgComp="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImgComp, plotlayoutImgComp, _, _=loadImgPlot(imgIntComp)
            btnOpticalDisable=false
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgComp=[traceImg]
            msgimgComp=""
        end
    end
    
    @onbutton imgPlusComp begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns())
        # Update the array of images listed in the public folder
        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        
        new_msi=increment_image(current_msiComp, msi_bmp)
        new_col_msi=increment_image(current_col_msiComp, col_msi_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_msiComp=new_msi
            current_col_msiComp=new_col_msi
            imgIntComp="/$(current_msiComp)?t=$(timestamp)"
            colorbarComp="/$(current_col_msiComp)?t=$(timestamp)"
    
            text_nmass=replace(current_msiComp, "MSI_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgimgComp="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImgComp, plotlayoutImgComp, _, _=loadImgPlot(imgIntComp)
            btnOpticalDisable=false
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgComp=[traceImg]
            msgimgComp=""
        end
    end
    
    @onbutton imgMinusTComp begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images with TrIQ filter listed in the public folder
        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    
        new_msi=decrement_image(current_triqComp, triq_bmp)
        new_col_msi=decrement_image(current_col_triqComp, col_triq_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_triqComp=new_msi
            current_col_triqComp=new_col_msi
            imgIntTComp="/$(current_triqComp)?t=$(timestamp)"
            colorbarTComp="/$(current_col_triqComp)?t=$(timestamp)"
    
            text_nmass=replace(current_triqComp, "TrIQ_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgtriqComp="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImgTComp, plotlayoutImgTComp, _, _=loadImgPlot(imgIntTComp)
            btnOpticalDisable=false
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgTComp=[traceImg]
            msgtriqComp=""
        end
    end
    
    @onbutton imgPlusTComp begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images with TrIQ filter listed in the public folder
        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    
        new_msi=increment_image(current_triqComp, triq_bmp)
        new_col_msi=increment_image(current_col_triqComp, col_triq_png)
        if new_msi!=nothing || new_col_msi!=nothing
            current_triqComp=new_msi
            current_col_triqComp=new_col_msi
            imgIntTComp="/$(current_triqComp)?t=$(timestamp)"
            colorbarTComp="/$(current_col_triqComp)?t=$(timestamp)"
    
            text_nmass=replace(current_triqComp, "TrIQ_" => "")
            text_nmass=replace(text_nmass, ".bmp" => "")
            msgtriqComp="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
            # Process the image in the function
            plotdataImgTComp, plotlayoutImgTComp, _, _=loadImgPlot(imgIntTComp)
            btnOpticalDisable=false
        else
            traceImg=PlotlyBase.heatmap(x=[], y=[])
            plotdataImgTComp=[traceImg]
            msgtriqComp=""
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
    end

    # Event detection for clicking on the spectrum plot
    @onchange data_click begin
        if selectedTab == "tab2"
            if !isempty(xSpectraMz)
                spectracoords=reshape(plotdata, 1, length(plotdata))
                # Extract x and y values from data_click 
                cursor_data=data_click["cursor"] 
                x_value=cursor_data["x"]
                y_value=cursor_data["y"]

                # Find the minimum x-value in spectracoords
                min_x_value=minimum([minimum(val[:x]) for val in spectracoords if !isempty(val[:x])])
                # Adjust x_value and spectracoords x-values to start from 0
                adjusted_x_value=x_value - min_x_value

                closest_distance=Inf
                for val in spectracoords
                    adjusted_x=val[:x] .- min_x_value
                    start_idx=findfirst(x -> x >= adjusted_x_value - 20, adjusted_x)
                    end_idx=findlast(x -> x <= adjusted_x_value + 20, adjusted_x)
                    
                    if start_idx !== nothing && end_idx !== nothing
                        for i in start_idx:end_idx
                            spectra_x=adjusted_x[i]
                            spectra_y=val[:y][i]
                            distance=sqrt((spectra_x - adjusted_x_value)^2 + (spectra_y - y_value)^2)
                            if distance < closest_distance
                                closest_distance=distance
                                Nmass=round(spectra_x + min_x_value, digits=2)
                            end
                        end
                    end
                end
                traceSpectra=PlotlyBase.scatter(x=xSpectraMz, y=ySpectraMz, mode="lines",name="Spectra",showlegend=false)
                trace2=PlotlyBase.scatter(x=[Nmass, Nmass],y=[0, maximum(ySpectraMz)],mode="lines",line=attr(color="red", width=0.5),name="<i>m/z</i> selected",showlegend=false)
                plotdata=[traceSpectra,trace2] # We add the data from spectra and the red line to the plot
            end
        elseif selectedTab == "tab1"
            cursor_data=data_click["cursor"] 
            xCoord=Int32(round(cursor_data["x"]))
            yCoord=Int32(round(cursor_data["y"]))
            if xCoord < 1
                xCoord=1
            elseif xCoord > imgWidth
                xCoord=imgWidth
            end
            if yCoord > -1
                yCoord=-1
            elseif yCoord < -imgHeight
                yCoord=-imgHeight
            end # Get the x and y values from the click of the cursor and make sure they don't exceed image proportions
            plotdataImgT=filter(trace -> !(get(trace, :name, "") in ["Line X", "Line Y"]), plotdataImgT)
            trace1, trace2=crossLinesPlot(xCoord, yCoord, imgWidth, -imgHeight)
            plotdataImgT=append!(plotdataImgT, [trace1, trace2])
        elseif selectedTab == "tab0"
            cursor_data=data_click["cursor"] 
            xCoord=Int32(round(cursor_data["x"]))
            yCoord=Int32(round(cursor_data["y"]))
            if xCoord < 1
                xCoord=1
            elseif xCoord > imgWidth
                xCoord=imgWidth
            end
            if yCoord > -1
                yCoord=-1
            elseif yCoord < -imgHeight
                yCoord=-imgHeight
            end # Get the x and y values from the click of the cursor and make sure they don't exceed image proportions
            plotdataImg=filter(trace -> !(get(trace, :name, "") in ["Line X", "Line Y","Optical"]), plotdataImg)
            trace1, trace2=crossLinesPlot(xCoord, yCoord, imgWidth, -imgHeight)
            plotdataImg=append!(plotdataImg, [trace1, trace2])
        end
    end

    @onbutton btnOptical begin
        imgRoute=pick_file(; filterlist="png,bmp,jpg,jpeg")
        selectedTab="tab0"
        plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
        if isnothing(imgRoute)
            msg="No optical image selected"
        else
            img=load(imgRoute)
            save("./public/css/imgOver.png",img)
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt,"/css/imgOver.png",imgTrans)
        end
        
    end

    @onbutton btnOpticalT begin
        imgRoute=pick_file(; filterlist="png,bmp,jpg,jpeg")
        selectedTab="tab1"
        plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
        if isnothing(imgRoute)
            msg="No optical image selected"
        else
            img=load(imgRoute)
            save("./public/css/imgOver.png",img)
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT,"/css/imgOver.png",imgTrans)
            opticalOverTriq=true
        end
    end

    @onchange imgTrans begin
        if !opticalOverTriq
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt,"/css/imgOver.png",imgTrans)
        else
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT,"/css/imgOver.png",imgTrans)
        end
    end

    @onchange opticalOverTriq begin
        if !opticalOverTriq
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt,"/css/imgOver.png",imgTrans)
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
            selectedTab="tab0"
        else
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT,"/css/imgOver.png",imgTrans)
            selectedTab="tab1"
        end
    end

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
