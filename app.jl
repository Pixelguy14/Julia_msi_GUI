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
@genietools

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded
rgb_ViridisPalette =reinterpret(ColorTypes.RGB24, ViridisPalette)

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
    @out btnStartDisable = false
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
    @in Main_Process = false # To generate images
    @in Generate_Plot = false # To generate plot
    @in progress = false
    @in progressPlot = false
    @in triqEnabled = false
    @in ImgPlus = false
    @in ImgMinus = false
    @in ImgPlusT = false
    @in ImgMinusT = false

    @out indeximg = 0
    @out indeximgTriq = 0

    @out lastimg = 0
    @out lastimgTriq = 0

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

    # WIP
    @out int_nmass = 0
    @out dec_nmass = 0.0

    # Interface Plot 
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

    traceSpectra = PlotlyBase.scatter(x=[], y=[], mode="lines")
    @out plotdata = [traceSpectra]
    @out plotlayout = layoutSpectra

    # == Reactive handlers ==
    # Reactive handlers watch a variable and execute a block of code when its value changes
    # The onbutton handler will set the variable to false after the block is executed
    """
    @onchange triqEnabled begin
        if !triqEnabled
            #triqProb = 0.0
            #triqColor = 0
        end
    end
    """

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
    
    @onbutton Main_Process begin
        progress = true # Start progress button animation
        btnStartDisable = true # We disable the button to avoid multiple requests
        indeximg = floor(Int, Nmass)
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
                        SaveBitmap(joinpath("public", "TrIQ_$(floor(Int, Nmass)).bmp"),TrIQ(slice, Int(triqColor), triqProb),ViridisPalette)
                        # Use timestamp to refresh image interface container
                        imgIntT = "/TrIQ_$(floor(Int, Nmass)).bmp?t=$(timestamp)"
                        msgtriq = "TrIQ image with the Nmass of $(floor(Int, Nmass))"
                        ticks = round.(range(0, stop = maximum(TrIQ(slice, Int(triqColor), triqProb)), length = 10), digits = 2)
                        Colorbar(fig[1, 1], colormap = rgb_ViridisPalette, limits = (0, maximum(TrIQ(slice, Int(triqColor), triqProb))),ticks = ticks, label = "Intensity")
                        save("public/colorbar_TrIQ_$(floor(Int, Nmass)).png", fig)
                        msg = "The file has been created successfully inside the 'public' folder of the app."
                        colorbarT = "/colorbar_TrIQ_$(floor(Int, Nmass)).png?t=$(timestamp)"
                    end
                else # If we don't use TrIQ
                    SaveBitmap(joinpath("public", "MSI_$(floor(Int, Nmass)).bmp"),IntQuant(slice),ViridisPalette)
                    # Use timestamp to refresh image interface container
                    imgInt = "/MSI_$(floor(Int, Nmass)).bmp?t=$(timestamp)"
                    msgimg = "image with the Nmass of $(floor(Int, Nmass))"
                    ticks = round.(range(0, stop = maximum(slice), length = 10), digits = 2)
                    Colorbar(fig[1, 1], colormap = rgb_ViridisPalette, limits = (0, maximum(slice)),ticks = ticks, label = "Intensity")
                    save("public/colorbar_MSI_$(floor(Int, Nmass)).png", fig)
                    msg = "The file has been created successfully inside the 'public' folder of the app."
                    colorbar = "/colorbar_MSI_$(floor(Int, Nmass)).png?t=$(timestamp)"
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

    @onbutton Generate_Plot begin
        msg = ""
        full_route = joinpath( file_route, file_name )
        progressPlot = true 
        if isfile(full_route) # Check if the file exists
            btnPlotDisable = false
            btnStartDisable = false
            full_routeMz = split( full_route, "." )[1] * ".mzML" # Splitting the route from imzml to mzml so the plotting can work
            if isfile(full_routeMz) && (full_routeMz2 == "" || full_routeMz2 != full_routeMz) # Check if the mzml exists
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

    # WIP THESE NEED TO CHANGE, A BETTER BLIND SEARCH AND ALLOWANCE FOR DECIMAL IMAGES
    # Also needs to implement timestamp to force reload of image and update if needed.
    @onbutton ImgMinus begin
        indeximg-=1
        while !isfile("public/MSI_$(indeximg).bmp") && indeximg > 0
            indeximg -= 1
        end
        if(indeximg <= 0) # If it doesn't find a lower value image
            indeximg = lastimg
        end
        imgInt = "/MSI_$(indeximg).bmp"
        colorbar = "/colorbar_MSI_$(indeximg).png"
        msgimg = "image with the Nmass of $(indeximg)"
        lastimg = indeximg
    end
    @onbutton ImgPlus begin
        indeximg+=1
        while !isfile("public/MSI_$(indeximg).bmp") && indeximg < 2001
            indeximg += 1
        end
        if(indeximg >= 2001)
            indeximg = lastimg
        end
        imgInt = "/MSI_$(indeximg).bmp"
        colorbar = "/colorbar_MSI_$(indeximg).png"
        msgimg = "image with the Nmass of $(indeximg)"
        lastimg = indeximg
    end

    @onbutton ImgMinusT begin
        indeximgTriq-=1
        while !isfile("public/TrIQ_$(indeximgTriq).bmp") && indeximgTriq > 0
            indeximgTriq -= 1
        end
        if(indeximgTriq <= 0) # If it doesn't find a lower value image
            indeximgTriq = lastimgTriq
        end
        imgIntT = "/TrIQ_$(indeximgTriq).bmp"
        colorbarT = "/colorbar_TrIQ_$(indeximgTriq).png"
        msgtriq = "TrIQ image with the Nmass of $(indeximgTriq)"
        lastimgTriq = indeximgTriq
    end
    @onbutton ImgPlusT begin
        indeximgTriq+=1
        while !isfile("public/TrIQ_$(indeximgTriq).bmp") && indeximgTriq < 2001
            indeximgTriq += 1
        end
        if(indeximgTriq >= 2001)
            indeximgTriq = lastimgTriq
        end
        imgIntT = "/TrIQ_$(indeximgTriq).bmp"
        colorbarT = "/colorbar_TrIQ_$(indeximgTriq).png"
        msgtriq = "TrIQ image with the Nmass of $(indeximgTriq)"
        lastimgTriq = indeximgTriq
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

