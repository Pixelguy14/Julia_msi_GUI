module App
# == Packages ==
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
using NativeFileDialog
@genietools

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded
rgb_ViridisPalette=reinterpret(ColorTypes.RGB24, ViridisPalette)

# == Search functions ==
function increment_image(current_image, image_list)
    if isempty(image_list)
        return nothing
    end
    current_index=findfirst(isequal(current_image), image_list)
    if current_index == nothing || current_index == length(image_list) || current_image === ""
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
    if current_index == nothing || current_index == 1 || current_image === ""
        return image_list[1]  # Return the current image if it's the first one or not found
    else
        return image_list[current_index - 1]  # Move to the previous image
    end
end

# == Reactive code ==
# reactive code to make the UI interactive
@app begin
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block
    #@out test="/test.bmp" #slash means it's getting the info from 'public' folder

    # Interface non Variables
    @out warning_fr=""
    @out btnStartDisable=true
    @out btnPlotDisable=false
    @in warning_msg=false
    @in CompareDialog=false

    # Interface Variables
    @in file_route=""
    @in file_name=""
    @in Nmass=0.0
    @in Tol=0.0
    @in triqProb=0.98
    @in triqColor=256

    # Interface Buttons
    @in btnSearch=false # To search for files in your device
    @in mainProcess=false # To generate images
    @in CompareBtn=false # to open dialog
    @in createSumPlot=false # To generate sum spectrum plot
    @in image3dPlot=false # To generate 3d plot based on current image
    @in triq3dPlot=false # To generate 3d plot based on current triq image
    @in imageCPlot=false # To generate contour plots of current image
    @in triqCPlot=false # To generate contour plots of current triq image
    @in progress=false
    @in progressPlot=false
    @in triqEnabled=false
    @in ImgPlus=false
    @in ImgMinus=false
    @in ImgPlusT=false
    @in ImgMinusT=false

    #TABS
    @out tabIDs=["tab0","tab1","tab2","tab3","tab4"]
    @out tabLabels=["Image", "TrIQ", "Spectrum Plot", "Topology Plot","Surface Plot"]
    @in selectedTab= "tab0"
    @out CompTabIDs=["tab0","tab1","tab2","tab3","tab4"]
    @out CompTabLabels=["Image", "TrIQ", "Spectrum Plot", "Topology Plot","Surface Plot"]
    @in CompSelectedTab= "tab0"

    # Interface Images
    @out imgInt="/.bmp" # image Interface
    @out imgIntT="/.bmp" # image Interface TrIQ
    @out colorbar="/.png"
    @out colorbarT="/.png"

    @out msg=""
    @out msgimg=""
    @out msgtriq=""

    @out full_route=""
    @out full_routeMz=""
    @out full_routeMz2=""

    # For the creation of images with a more specific mass charge
    @out text_nmass=""

    # For image search
    # Image lists we apply a filter that searches specific type of images into our public folder, then we sort it in a "numerical" order
    @in msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    @in triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    # Set current image for the list
    @out current_msi=""
    @out current_col_msi=""
    @out current_triq=""
    @out current_col_triq=""

    @out indeximg=0
    @out indeximgTriq=0

    @out lastimg=0
    @out lastimgTriq=0

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
    # Dummy 2D surface plot
    traceSpectra=PlotlyBase.scatter(x=[], y=[], mode="lines")
    # Create conection to frontend
    @out plotdata=[traceSpectra]
    @out plotlayout=layoutSpectra

    # Interface Plot Surface
    layoutContour=PlotlyBase.Layout(
        title="2D Topographic Map",
        xaxis=PlotlyBase.attr(
            title="X",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Y",
            showgrid=true
        )
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
    
    # println("3D trace defined: ", trace3D)

    # == Reactive handlers ==
    # Reactive handlers watch a variable and execute a block of code when its value changes
    # The onbutton handler will set the variable to false after the block is executed

    @onbutton btnSearch begin
        full_route=pick_file(; filterlist="imzML")
        if isnothing(full_route)
            #println("No file selected")
            msg="No file selected"
            warning_msg=true
            btnStartDisable=true
        else
            #println("Selected file path: ", full_route)
            btnStartDisable=false
            btnPlotDisable=false
            msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
            col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
            triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
            col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        end
    end

    """
    @onchange Nmass begin
        indeximg=floor(Int, Nmass)
        indeximgTriq=floor(Int, Nmass)
        lastimg=floor(Int, Nmass)
        lastimgTriq=floor(Int, Nmass)
    end
    """
    
    @onbutton mainProcess begin
        progress=true # Start progress button animation
        btnStartDisable=true # We disable the button to avoid multiple requests
        btnPlotDisable=true
        indeximg=floor(Int, Nmass)
        text_nmass=replace(string(Nmass), "." => "_")
        #full_route=joinpath(file_route, file_name)
        if isfile(full_route) && Nmass > 0 && Tol > 0 && Tol <= 1
            msg="File exists, Nmass=$(Nmass) Tol=$(Tol). Loading file will begin, please be patient."
            try
                spectra=LoadImzml(full_route)
                msg="File loaded. Creating Spectra with the specific mass and tolerance, please be patient."
                slice=GetSlice(spectra, Nmass, Tol)
                fig=CairoMakie.Figure(size=(100, 200)) # Container
                # Append a query string to force the image to refresh 
                timestamp=string(time_ns()) 
                if triqEnabled # If we have TrIQ
                    if triqColor < 1 || triqColor > 256 ||triqProb < 0 || triqProb > 1
                        msg="Incorrect TrIQ values, please adjust accordingly and try again."
                        warning_msg=true
                    else
                        image_path=joinpath("./public", "TrIQ_$(text_nmass).bmp")
                        SaveBitmap(joinpath("public", "TrIQ_$(text_nmass).bmp"),TrIQ(slice, Int(triqColor), triqProb),ViridisPalette)
                        # Flip te image vertically then save it again
                        img=load(image_path)
                        if size(img, 1) > size(img, 2) # fix to taller images 
                            img=reverse(permutedims(img, (2, 1)), dims=1)
                        end
                        flipped_img=reverse(img, dims=1)
                        save(image_path, flipped_img)
                        # Use timestamp to refresh image interface container
                        imgIntT="/TrIQ_$(text_nmass).bmp?t=$(timestamp)"
                        # Get current image 
                        current_triq="TrIQ_$(text_nmass).bmp"
                        msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
                        # Create colorbar 
                        ticks=round.(range(0, stop=maximum(TrIQ(slice, Int(triqColor), triqProb)), length=10), digits=2)
                        Colorbar(fig[1, 1], colormap=rgb_ViridisPalette, limits=(0, maximum(TrIQ(slice, Int(triqColor), triqProb))),ticks=ticks, label="Intensity")
                        save("public/colorbar_TrIQ_$(text_nmass).png", fig)
                        colorbarT="/colorbar_TrIQ_$(text_nmass).png?t=$(timestamp)"
                        # Get current colorbar 
                        current_col_triq="colorbar_TrIQ_$(text_nmass).png"
                        # We update the directory to include the new placed images.
                        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
                        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
                        msg="The file has been created successfully inside the 'public' folder of the app."
                        selectedTab="tab1"
                        #println("all msi in folder=",triq_bmp)
                        #println("all col msi in folder= ",col_triq_png)
                    end
                else # If we don't use TrIQ
                    image_path=joinpath("./public", "MSI_$(text_nmass).bmp")
                    SaveBitmap(joinpath("public", "MSI_$(text_nmass).bmp"),IntQuant(slice),ViridisPalette)
                    # Flip te image vertically then save it again
                    img=load(image_path)
                    if size(img, 1) > size(img, 2) # fix to taller images 
                        img=reverse(permutedims(img, (2, 1)), dims=1) 
                    end
                    flipped_img=reverse(img, dims=1)
                    save(image_path, flipped_img)
                    # Use timestamp to refresh image interface container
                    imgInt="/MSI_$(text_nmass).bmp?t=$(timestamp)"
                    # Get current image 
                    current_msi="MSI_$(text_nmass).bmp"
                    msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
                    # Create colorbar 
                    ticks=round.(range(0, stop=maximum(slice), length=10), digits=2)
                    Colorbar(fig[1, 1], colormap=rgb_ViridisPalette, limits=(0, maximum(slice)),ticks=ticks, label="Intensity")
                    save("public/colorbar_MSI_$(text_nmass).png", fig)
                    colorbar="/colorbar_MSI_$(text_nmass).png?t=$(timestamp)"
                    # Get current colorbar 
                    current_col_msi="colorbar_MSI_$(text_nmass).png"
                    # We update the directory to include the new placed images.
                    msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
                    col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
                    selectedTab="tab0"
                    msg="The file has been created successfully inside the 'public' folder of the app."
                    #println("all msi in folder=",msi_bmp)
                    #println("all col msi in folder= ",col_msi_png)
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
        progress=false
    end

    @onbutton createSumPlot begin
        msg="Sum spectrum plot selected"
        #full_route=joinpath( file_route, file_name )
        if isfile(full_route) # Check if the file exists
            btnPlotDisable=false
            btnStartDisable=false
            full_routeMz=split( full_route, "." )[1] * ".mzML" # Splitting the route from imzml to mzml so the plotting can work
            if isfile(full_routeMz) && (full_routeMz2 == "" || full_routeMz2 != full_routeMz) # Check if the mzml exists
                println("I'm working as intended")
                progressPlot=true
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
                    )
                )
                # dims=size(spectraMz)
                # scansMax=dims[2] # we get the total of scansMax
                # traceSpectra=PlotlyBase.scatter(x=spectraMz[1, 1], y=spectraMz[2, 1], mode="lines")
                traceSpectra=PlotlyBase.scatter(x=mean(spectraMz[1,:]), y=mean(spectraMz[2,:]), mode="lines")
                plotdata=[traceSpectra] # We add the data from spectra to the plot
                plotlayout=layoutSpectra
                spectraMz=nothing # Important for memory cleaning
                GC.gc() # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                end
                selectedTab="tab2"
                msg="Plot loaded."
                full_routeMz2=full_routeMz # To avoid creating the plot if its the same file read as before
            else
                msg="the mzML file was not found"
                warning_msg=true
            end
        else
            msg="is not an imzML file"
            warning_msg=true
        end
        progressPlot=false
        btnPlotDisable=false
        btnStartDisable=false
    end

    # Image loaders based on the position of the current image (increment and decrement for both normal and filter)
    # And a pre-generated list from all image files from /public folder
    @onbutton ImgMinus begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        # Update the array of images listed in the public folder
            msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
            col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        
        new_msi=decrement_image(current_msi, msi_bmp)
        new_col_msi=decrement_image(current_col_msi, col_msi_png)
        
        current_msi=new_msi
        current_col_msi=new_col_msi
        imgInt="/$(current_msi)?t=$(timestamp)"
        colorbar="/$(current_col_msi)?t=$(timestamp)"

        text_nmass=replace(current_msi, "MSI_" => "")
        text_nmass=replace(text_nmass, ".bmp" => "")
        msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
    end
    @onbutton ImgPlus begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns())
        # Update the array of images listed in the public folder
            msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
            col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        
        new_msi=increment_image(current_msi, msi_bmp)
        new_col_msi=increment_image(current_col_msi, col_msi_png)

        current_msi=new_msi
        current_col_msi=new_col_msi
        imgInt="/$(current_msi)?t=$(timestamp)"
        colorbar="/$(current_col_msi)?t=$(timestamp)"

        text_nmass=replace(current_msi, "MSI_" => "")
        text_nmass=replace(text_nmass, ".bmp" => "")
        msgimg="Image with the Nmass of $(replace(text_nmass, "_" => "."))"
    end

    @onbutton ImgMinusT begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        new_msi=decrement_image(current_triq, triq_bmp)
        new_col_msi=decrement_image(current_col_triq, col_triq_png)
        # Update the array of images with TrIQ filter listed in the public folder
        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)

        current_triq=new_msi
        current_col_triq=new_col_msi
        imgIntT="/$(current_triq)?t=$(timestamp)"
        colorbarT="/$(current_col_triq)?t=$(timestamp)"

        text_nmass=replace(current_triq, "TrIQ_" => "")
        text_nmass=replace(text_nmass, ".bmp" => "")
        msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
        
    end
    @onbutton ImgPlusT begin
        # Append a query string to force the image to refresh 
        timestamp=string(time_ns()) 
        new_msi=increment_image(current_triq, triq_bmp)
        new_col_msi=increment_image(current_col_triq, col_triq_png)
        # Update the array of images with TrIQ filter listed in the public folder
        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        
        current_triq=new_msi
        current_col_triq=new_col_msi
        imgIntT="/$(current_triq)"
        colorbarT="/$(current_col_triq)"

        text_nmass=replace(current_triq, "TrIQ_" => "")
        text_nmass=replace(text_nmass, ".bmp" => "")
        msgtriq="TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
    end

    # 3d plot
    @onbutton image3dPlot begin
        msg="Image 3D plot selected"
        cleaned_imgInt=replace(imgInt, r"\?.*" => "")
        cleaned_imgInt=lstrip(cleaned_imgInt, '/')
        var=joinpath( "./public", cleaned_imgInt )

        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            try
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
                if size(elevation_smoothed, 1) < size(elevation_smoothed, 2)
                    # Transpose the elevation_smoothed array if Y axis is longer than X axis to fix chopping
                    elevation_smoothed=transpose(elevation_smoothed)
                    Y=-Y
                end

                trace3D=PlotlyBase.surface(x=X[1, :], y=Y[:, 1], z=elevation_smoothed,
                contours_z=attr(
                    show=true,
                    usecolormap=true,
                    highlightcolor="limegreen",
                    project_z=true
                ), colorscale="Viridis")
                plotdata3d=[trace3D] # We add the data from the image to the plot
                plotlayout3d=layout3D # we update the style of the plot to fit the image.
                spectraMz=nothing # Important for memory cleaning
                GC.gc() # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                end
                selectedTab="tab4"
                msg="Plot loaded."
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
    end
    # 3d plot for TrIQ
    @onbutton triq3dPlot begin
        msg="TrIQ 3D plot selected"
        cleaned_imgIntT=replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT=lstrip(cleaned_imgIntT, '/')
        var=joinpath( "./public", cleaned_imgIntT )

        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            try
                img=load(var)
                img_gray=Gray.(img) # Convert to grayscale
                img_array=Array(img_gray)
                elevation=Float32.(Array(img_gray)) ./ 255.0 # Normalize between 0 and 1
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
                    title="3D Surface Plot",
                    scene=attr(
                        xaxis_nticks=x_nticks,
                        yaxis_nticks=y_nticks,
                        zaxis_nticks=z_nticks, 
                        camera=attr(eye=attr(x=0, y=-1, z=0.5)), 
                        aspectratio=aspect_ratio
                    )
                )
                if size(elevation_smoothed, 1) < size(elevation_smoothed, 2)
                    # Transpose the elevation_smoothed array if Y axis is longer than X axis to fix chopping
                    elevation_smoothed=transpose(elevation_smoothed)
                    Y=-Y
                end

                trace3D=PlotlyBase.surface(x=X[1, :], y=Y[:, 1], z=elevation_smoothed,
                contours_z=attr(
                    show=true,
                    usecolormap=true,
                    highlightcolor="limegreen",
                    project_z=true
                ), colorscale="Viridis")
                plotdata3d=[trace3D] # We add the data from the image to the plot
                plotlayout3d=layout3D # we update the style of the plot to fit the image.
                spectraMz=nothing # Important for memory cleaning
                GC.gc() # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                end
                selectedTab="tab4"
                msg="Plot loaded."
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
    end

    # Contour 2d plot
    @onbutton imageCPlot begin
        msg="Image 2D plot selected"
        cleaned_imgInt=replace(imgInt, r"\?.*" => "")
        cleaned_imgInt=lstrip(cleaned_imgInt, '/')
        var=joinpath("./public", cleaned_imgInt)
    
        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            try
                img=load(var)
                # Convert to grayscale
                img_gray=Gray.(img)
                img_array=Array(img_gray)
                elevation=Float32.(Array(img_gray)) ./ 255.0  # Normalize between 0 and 1
    
                # Smooth the image
                sigma=3.0
                kernel=Kernel.gaussian(sigma)
                elevation_smoothed=imfilter(elevation, kernel)
    
                # Create the X, Y meshgrid coordinates
                x=1:size(elevation_smoothed, 2)
                y=1:size(elevation_smoothed, 1)
                X=repeat(reshape(x, 1, length(x)), length(y), 1)
                Y=repeat(reshape(y, length(y), 1), 1, length(x))
                
                layoutContour=PlotlyBase.Layout(
                    title="2D Topographic Map",
                    xaxis_title="X",
                    yaxis_title="Y"
                    
                )
                traceContour=PlotlyBase.contour(
                    z=elevation_smoothed,
                    x=X[1, :],  # Use the first row
                    y=-Y[:, 1],  # Use the first column
                    contours_coloring="Viridis",
                    colorscale="Viridis"
                )
                plotdataC=[traceContour]
                plotlayoutC=layoutContour
    
                spectraMz=nothing  # Important for memory cleaning
                GC.gc()  # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0)  # Ensure Julia returns the freed memory to OS
                end
                selectedTab="tab3"
                msg="Plot loaded."
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
    end
    # Contour 2d plot for TrIQ
    @onbutton triqCPlot begin
        msg="Image 2D plot selected"
        cleaned_imgIntT=replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT=lstrip(cleaned_imgIntT, '/')
        var=joinpath("./public", cleaned_imgIntT)
    
        if isfile(var)
            progressPlot=true
            btnPlotDisable=true
            btnStartDisable=true
            try
                img=load(var)
                # Convert to grayscale
                img_gray=Gray.(img)
                img_array=Array(img_gray)
                elevation=Float32.(Array(img_gray)) ./ 255.0  # Normalize between 0 and 1
    
                # Smooth the image
                sigma=3.0
                kernel=Kernel.gaussian(sigma)
                elevation_smoothed=imfilter(elevation, kernel)
    
                # Create the X, Y meshgrid coordinates
                x=1:size(elevation_smoothed, 2)
                y=1:size(elevation_smoothed, 1)
                X=repeat(reshape(x, 1, length(x)), length(y), 1)
                Y=repeat(reshape(y, length(y), 1), 1, length(x))
                
                layoutContour=PlotlyBase.Layout(
                    title="2D Topographic Map",
                    xaxis_title="X",
                    yaxis_title="Y",
                    
                )
                traceContour=PlotlyBase.contour(
                    z=elevation_smoothed,
                    x=X[1, :],  # Use the first row
                    y=-Y[:, 1],  # Use the first column
                    contours_coloring="Viridis",
                    colorscale="Viridis"
                )
                plotdataC=[traceContour]
                plotlayoutC=layoutContour
    
                spectraMz=nothing  # Important for memory cleaning
                GC.gc()  # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0)  # Ensure Julia returns the freed memory to OS
                end
                selectedTab="tab3"
                msg="Plot loaded."
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
    end

    @onbutton CompareBtn begin
        CompareDialog=true
    end

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

