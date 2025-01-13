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
# using ImageMagick
@genietools

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded
rgb_ViridisPalette = reinterpret(ColorTypes.RGB24, ViridisPalette)

# == Search functions ==
function increment_image(current_image, image_list)
    current_index = findfirst(isequal(current_image), image_list)
    if current_index == nothing || current_index == length(image_list) || current_image === ""
        return image_list[length(image_list)]  # Return the current image if it's the last one or not found
    else
        return image_list[current_index + 1]  # Move to the next image
    end
end

function decrement_image(current_image, image_list)
    current_index = findfirst(isequal(current_image), image_list)
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
    #@out test = "/test.bmp" #slash means it's getting the info from 'public' folder

    # Interface non Variables
    @out warning_fr = ""
    @out btnStartDisable = true
    @out btnPlotDisable = true
    @in warning_msg = false

    # Interface Variables
    @in file_route = ""
    @in file_name = ""
    @in Nmass = 0.0
    @in Tol = 0.0
    @in triqProb = 0.0
    @in triqColor = 0

    # Interface Buttons
    @in mainProcess = false # To generate images
    @in createSumPlot = false # To generate sum spectrum plot
    @in image3dPlot = false # To generate 3d plot based on current image
    @in triq3dPlot = false # To generate 3d plot based on current triq image
    @in progress = false
    @in progressPlot = false
    @in triqEnabled = false
    @in ImgPlus = false
    @in ImgMinus = false
    @in ImgPlusT = false
    @in ImgMinusT = false

    # Interface Images
    @out imgInt = "/.bmp" # image Interface
    @out imgIntT = "/.bmp" # image Interface TrIQ
    @out colorbar = "/.png"
    @out colorbarT = "/.png"

    @out msg = ""
    @out msgimg = ""
    @out msgtriq = ""

    @out full_route = ""
    @out full_routeMz = ""
    @out full_routeMz2 = ""

    # For the creation of images with a more specific mass charge
    @out text_nmass = ""

    # For image search
    # Image lists we apply a filter that searches specific type of images into our public folder, then we sort it in a "numerical" order
    @in msi_bmp = sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_msi_png = sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    @in triq_bmp = sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_triq_png = sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    # Set current image for the list
    @out current_msi = ""
    @out current_col_msi = ""
    @out current_triq = ""
    @out current_col_triq = ""

    @out indeximg = 0
    @out indeximgTriq = 0

    @out lastimg = 0
    @out lastimgTriq = 0

    # Interface Plot 2d
    layoutSpectra = PlotlyBase.Layout(
        title = "SUM Spectrum plot",
        xaxis = PlotlyBase.attr(
            title = "<i>m/z</i>",
            showgrid = true
        ),
        yaxis = PlotlyBase.attr(
            title = "Intensity",
            showgrid = true
        )
    )
    # Dummy 2D surface plot
    traceSpectra = PlotlyBase.scatter(x=[], y=[], mode="lines")
    # Create conection to frontend
    @out plotdata = [traceSpectra]
    @out plotlayout = layoutSpectra

    # Interface Plot 3d
    # Define the layout for the 3D plot
    layout3D = PlotlyBase.Layout(
        title = "3D Surface Plot",
        scene = attr(
            xaxis_title = "X",
            yaxis_title = "Y",
            zaxis_title = "Z"
        )
    )

    # Dummy 3D surface plot
    x = 1:10
    y = 1:10
    z = [sin(i * j / 10) for i in x, j in y]
    trace3D = PlotlyBase.surface(x=x, y=y, z=z,
                                contours_z=attr(
                                    show=true,
                                    usecolormap=true,
                                    highlightcolor="limegreen",
                                    project_z=true
                                ), colorscale="Viridis")
    # Create conection to frontend
    @out plotdata3d = [trace3D]
    @out plotlayout3d = layout3D
    
    # println("3D trace defined: ", trace3D)

    # == Reactive handlers ==
    # Reactive handlers watch a variable and execute a block of code when its value changes
    # The onbutton handler will set the variable to false after the block is executed
    @onchange triqEnabled begin
        # just in case that they remove an image manually, we use the least intensive button to get all images and re-order them.
        msi_bmp = sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_msi_png = sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
        triq_bmp = sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
        col_triq_png = sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    end
    

    @onchange file_name begin
        msg = ""
        progress = false
        progressPlot = false
        try
            if contains(file_name,".imzML")
                warning_fr = ""
                full_route = joinpath( file_route, file_name )
                if isfile(full_route) # Check if the file exists
                    btnPlotDisable = false
                    btnStartDisable = false
                else
                    warning_fr = "is not an imzML file"
                end
            else
                btnPlotDisable = true
                btnStartDisable = true
                full_route = "/"
                warning_fr = "is not an imzML or mzML file"
            end
        catch e
            msg = "There was an error, please verify the file and try again. $(e)"
        end
    end
    """
    @onchange Nmass begin
        indeximg = floor(Int, Nmass)
        indeximgTriq = floor(Int, Nmass)
        lastimg = floor(Int, Nmass)
        lastimgTriq = floor(Int, Nmass)
    end
    """
    
    @onbutton mainProcess begin
        progress = true # Start progress button animation
        btnStartDisable = true # We disable the button to avoid multiple requests
        indeximg = floor(Int, Nmass)
        text_nmass = replace(string(Nmass), "." => "_")
        full_route = joinpath(file_route, file_name)
        if isfile(full_route) && Nmass > 0 && Tol > 0 && Tol <= 1
            msg = "File exists, Nmass=$(Nmass) Tol=$(Tol). Loading file will begin, please be patient."
            try
                spectra = LoadImzml(full_route)
                msg = "File loaded. Creating Spectra with the specific mass and tolerance, please be patient."
                slice = GetSlice(spectra, Nmass, Tol)
                fig = CairoMakie.Figure(size = (100, 200)) # Container
                # Append a query string to force the image to refresh 
                timestamp = string(time_ns()) 
                if triqEnabled # If we have TrIQ
                    if triqColor < 1 || triqColor > 256 ||triqProb < 0 || triqProb > 1
                        msg = "Incorrect TrIQ values, please adjust accordingly and try again."
                        warning_msg = true
                    else
                        SaveBitmap(joinpath("public", "TrIQ_$(text_nmass).bmp"),TrIQ(slice, Int(triqColor), triqProb),ViridisPalette)
                        # Use timestamp to refresh image interface container
                        imgIntT = "/TrIQ_$(text_nmass).bmp?t=$(timestamp)"
                        # Get current image 
                        current_triq = "TrIQ_$(text_nmass).bmp"
                        msgtriq = "TrIQ image with the Nmass of $(replace(text_nmass, "_" => "."))"
                        # Create colorbar 
                        ticks = round.(range(0, stop = maximum(TrIQ(slice, Int(triqColor), triqProb)), length = 10), digits = 2)
                        Colorbar(fig[1, 1], colormap = rgb_ViridisPalette, limits = (0, maximum(TrIQ(slice, Int(triqColor), triqProb))),ticks = ticks, label = "Intensity")
                        save("public/colorbar_TrIQ_$(text_nmass).png", fig)
                        colorbarT = "/colorbar_TrIQ_$(text_nmass).png?t=$(timestamp)"
                        # Get current colorbar 
                        current_col_triq = "colorbar_TrIQ_$(text_nmass).png"
                        # We update the directory to include the new placed images.
                        triq_bmp = sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
                        col_triq_png = sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)

                        msg = "The file has been created successfully inside the 'public' folder of the app."
                        #println("all msi in folder = ",triq_bmp)
                        #println("all col msi in folder= ",col_triq_png)
                    end
                else # If we don't use TrIQ
                    SaveBitmap(joinpath("public", "MSI_$(text_nmass).bmp"),IntQuant(slice),ViridisPalette)
                    # Use timestamp to refresh image interface container
                    imgInt = "/MSI_$(text_nmass).bmp?t=$(timestamp)"
                    # Get current image 
                    current_msi = "MSI_$(text_nmass).bmp"
                    msgimg = "image with the Nmass of $(replace(text_nmass, "_" => "."))"
                    # Create colorbar 
                    ticks = round.(range(0, stop = maximum(slice), length = 10), digits = 2)
                    Colorbar(fig[1, 1], colormap = rgb_ViridisPalette, limits = (0, maximum(slice)),ticks = ticks, label = "Intensity")
                    save("public/colorbar_MSI_$(text_nmass).png", fig)
                    colorbar = "/colorbar_MSI_$(text_nmass).png?t=$(timestamp)"
                    # Get current colorbar 
                    current_col_msi = "colorbar_MSI_$(text_nmass).png"
                    # We update the directory to include the new placed images.
                    msi_bmp = sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
                    col_msi_png = sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)

                    msg = "The file has been created successfully inside the 'public' folder of the app."
                    #println("all msi in folder = ",msi_bmp)
                    #println("all col msi in folder= ",col_msi_png)
                end
            catch e
                msg = "There was an error loading the ImzML file, please verify the file accordingly and try again. $(e)"
                warning_msg = true
            end
        else
            msg = "File does not exist or a parameter is incorrect, please try again."
            warning_msg = true
        end
        spectra = nothing # Important for memory cleaning
        slice = nothing
        GC.gc() # Trigger garbage collection
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
        end
        btnStartDisable = false
        progress = false
    end

    @onbutton createSumPlot begin
        msg = "Sum spectrum plot selected"
        full_route = joinpath( file_route, file_name )
        progressPlot = true
        if isfile(full_route) # Check if the file exists
            btnPlotDisable = false
            btnStartDisable = false
            full_routeMz = split( full_route, "." )[1] * ".mzML" # Splitting the route from imzml to mzml so the plotting can work
            if isfile(full_routeMz) && (full_routeMz2 == "" || full_routeMz2 != full_routeMz) # Check if the mzml exists
                println("I'm working as intended")
                btnPlotDisable = true
                msg = "Loading plot..."
                spectraMz = LoadMzml(full_routeMz)
                # dims = size(spectraMz)
                # scansMax = dims[2] # we get the total of scansMax
                # traceSpectra = PlotlyBase.scatter(x = spectraMz[1, 1], y = spectraMz[2, 1], mode="lines")
                traceSpectra = PlotlyBase.scatter(x = mean(spectraMz[1,:]), y = mean(spectraMz[2,:]), mode="lines")
                plotdata = [traceSpectra] # We add the data from spectra to the plot
                spectraMz = nothing # Important for memory cleaning
                GC.gc() # Trigger garbage collection
                if Sys.islinux()
                    ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                end
                msg = "Plot loaded."
                btnPlotDisable = false
                full_routeMz2 = full_routeMz # To avoid creating the plot if its the same file read as before
            end
        else
            msg = "is not an imzML file"
            warning_msg = true
        end
        progressPlot = false
    end

    # Image loaders based on the position of the current image (increment and decrement for both normal and filter)
    # And a pre-generated list from all image files from /public folder
    @onbutton ImgMinus begin
        # Append a query string to force the image to refresh 
        timestamp = string(time_ns()) 
        new_msi = decrement_image(current_msi, msi_bmp)
        new_col_msi = decrement_image(current_col_msi, col_msi_png)
        
        current_msi = new_msi
        current_col_msi = new_col_msi
        imgInt = "/$(current_msi)?t=$(timestamp)"
        colorbar = "/$(current_col_msi)?t=$(timestamp)"

        text_nmass = replace(current_msi, "MSI_" => "")
        text_nmass = replace(text_nmass, ".bmp" => "")
        msgimg = "image with the Nmass of $(replace(text_nmass, "_" => "."))"
    end
    @onbutton ImgPlus begin
        # Append a query string to force the image to refresh 
        timestamp = string(time_ns()) 
        new_msi = increment_image(current_msi, msi_bmp)
        new_col_msi = increment_image(current_col_msi, col_msi_png)

        current_msi = new_msi
        current_col_msi = new_col_msi
        imgInt = "/$(current_msi)?t=$(timestamp)"
        colorbar = "/$(current_col_msi)?t=$(timestamp)"

        text_nmass = replace(current_msi, "MSI_" => "")
        text_nmass = replace(text_nmass, ".bmp" => "")
        msgimg = "image with the Nmass of $(replace(text_nmass, "_" => "."))"
    end

    @onbutton ImgMinusT begin
        # Append a query string to force the image to refresh 
        timestamp = string(time_ns()) 
        new_msi = decrement_image(current_triq, triq_bmp)
        new_col_msi = decrement_image(current_col_triq, col_triq_png)

        current_triq = new_msi
        current_col_triq = new_col_msi
        imgIntT = "/$(current_triq)?t=$(timestamp)"
        colorbarT = "/$(current_col_triq)?t=$(timestamp)"

        text_nmass = replace(current_triq, "TrIQ_" => "")
        text_nmass = replace(text_nmass, ".bmp" => "")
        msgtriq = "image with the Nmass of $(replace(text_nmass, "_" => "."))"
        
    end
    @onbutton ImgPlusT begin
        # Append a query string to force the image to refresh 
        timestamp = string(time_ns()) 
        new_msi = increment_image(current_triq, triq_bmp)
        new_col_msi = increment_image(current_col_triq, col_triq_png)

        current_triq = new_msi
        current_col_triq = new_col_msi
        imgIntT = "/$(current_triq)"
        colorbarT = "/$(current_col_triq)"

        text_nmass = replace(current_triq, "TrIQ_" => "")
        text_nmass = replace(text_nmass, ".bmp" => "")
        msgtriq = "image with the Nmass of $(replace(text_nmass, "_" => "."))"
    end

    # 3d plot 
    @onbutton image3dPlot begin
        msg = "Image 3D plot selected"
        cleaned_imgInt = replace(imgInt, r"\?.*" => "")
        cleaned_imgInt = lstrip(cleaned_imgInt, '/')
        var1 = joinpath( "public", cleaned_imgInt )
        if isfile(var1)
            try
                img = FileIO.load(File{DataFormat{:BMP}}(var1))
                img_gray = Gray.(img) # Convert to grayscale 
                img_array = Array(img_gray)
            catch e 
                msg = "Failed to load and process image: $e" 
                warning_msg = true 
                println(msg) 
            end
            elevation = img_array / 255.0 # Normalize between 0 and 1
            # Smooth the image 
            sigma = 4.0 
            elevation_smoothed = imfilter(elevation, Kernel.gaussian(sigma)) 
            # Create the X, Y meshgrid coordinates 
            x = 1:size(elevation_smoothed, 2) 
            y = 1:size(elevation_smoothed, 1) 
            X, Y = ndgrid(x, y)
            trace3D = PlotlyBase.scatter(x = X[1, :], y =Y[:, 1], z=elevation_smoothed, colorscale="Viridis")
            plotdata3d = [trace3D] # We add the data from the image to the plot
            spectraMz = nothing # Important for memory cleaning
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end
             msg = "Plot loaded."
        else
            msg = "image could not be 3d plotted"
            warning_msg = true
        end
    end

    @onbutton triq3dPlot begin
        msg = "TrIQ 3D plot selected"
        cleaned_imgIntT = replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT = lstrip(cleaned_imgIntT, '/')
        var1 = joinpath( "public", cleaned_imgIntT )
        if isfile(var1)
            try
                img = FileIO.load(File{DataFormat{:BMP}}(var1))
                img_gray = Gray.(img) # Convert to grayscale 
                img_array = Array(img_gray)
            catch e 
                msg = "Failed to load and process image: $e" 
                warning_msg = true 
                println(msg) 
            end
            elevation = img_array / 255.0 # Normalize between 0 and 1
            # Smooth the image 
            sigma = 4.0 
            elevation_smoothed = imfilter(elevation, Kernel.gaussian(sigma)) 
            # Create the X, Y meshgrid coordinates 
            x = 1:size(elevation_smoothed, 2) 
            y = 1:size(elevation_smoothed, 1) 
            X, Y = ndgrid(x, y)
            trace3D = PlotlyBase.scatter(x = X[1, :], y =Y[:, 1], z=elevation_smoothed, colorscale="Viridis")
            plotdata3d = [trace3D] # We add the data from the image to the plot
            spectraMz = nothing # Important for memory cleaning
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end
             msg = "Plot loaded."
        else
            msg = "image could not be 3d plotted"
            warning_msg = true
        end
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
    @private table = DataFrame(a = 1:10, b = 10:19, c = 20:29)

=#

