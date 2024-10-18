module App
# == Packages ==
using GenieFramework # set up Genie development environment.
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
    @in file_route = ""
    @in file_name = ""
    @out warning_fr = ""
    @in Nmass = 0.0
    @in Tol = 0.0
    @in triqEnabled = false
    @out Disab_btn = false
    @in triqProb = 0.0
    @in triqColor = 0
    @in Main_Process = false
    @out indeximg = 0
    @out indeximgTriq = 0
    @out lastimg = 0
    @out lastimgTriq = 0
    @in ImgPlus = false
    @in ImgMinus = false
    @out test = "/.bmp"
    @in ImgPlusT = false
    @in ImgMinusT = false
    @out testT = "/.bmp"
    @out msg = ""
    @out msgimg = ""
    @out msgtriq = ""
    @out full_route = ""
    @out full_routeMz = ""
    @out full_routeMz2 = ""
    @out colorbar = "/.png"
    @out colorbarT = "/.png"
    layoutSpectra = PlotlyBase.Layout(
        title = "SUM Spectrum plot",
        xaxis = PlotlyBase.attr(
            title = "<i>m/z</i>",
            showgrid = true
        ),
        yaxis = PlotlyBase.attr(
            title = "Intensity",
            showgrid = true
        ),
        width = 900,
        height = 500
        )
    traceSpectra = PlotlyBase.scatter(x=[], y=[], mode="lines")
    @out plotdata = [traceSpectra]
    @out plotlayout = layoutSpectra

    # == Reactive handlers ==
    # reactive handlers watch a variable and execute a block of code when
    # its value changes
    @onchange triqEnabled begin
        if !triqEnabled
            #triqProb = 0.0
            #triqColor = 0
        end
    end
    @onchange file_name begin
        msg = ""
        try
            if contains(file_name,".imzML")
                warning_fr = ""
                full_route = joinpath( file_route, file_name )
                if isfile(full_route) # check if the file exists
                    full_routeMz = split( full_route, "." )[1] * ".mzML" # Splitting the route from imzml to mzml so the plotting can work
                    if isfile(full_routeMz) && (full_routeMz2 == "" || full_routeMz2 != full_routeMz) # check if there is an mzML file around
                        Disab_btn = true
                        warning_fr = "Loading plot..."
                        spectraMz = LoadMzml(full_routeMz)
                        # dims = size(spectraMz)
                        # scansMax = dims[2] # we get the total of scansMax
                        # traceSpectra = PlotlyBase.scatter(x = spectraMz[1, 1], y = spectraMz[2, 1], mode="lines")
                        traceSpectra = PlotlyBase.scatter(x = mean(spectraMz[1,:]), y = mean(spectraMz[2,:]), mode="lines")
                        plotdata = [traceSpectra] # we add the data of spectra to the plot
                        spectraMz = nothing # Important for memory cleaning
                        GC.gc() # Trigger garbage collection
                        if Sys.islinux()
                            ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                        end
                        warning_fr = "Plot loaded."
                        Disab_btn = false
                        full_routeMz2 = full_routeMz # to avoid creating the plot if its the same file read as before
                    end
                else
                    warning_fr = "is not an imzML file"
                end
            elseif contains(file_name,".mzML")
                full_routeMz = joinpath( file_route, file_name )
                warning_fr = "$(full_routeMz)"
                if isfile(full_routeMz) && (full_routeMz2 == "" || full_routeMz2 != full_routeMz) # check if there is an mzML file around # check if the file exists
                    Disab_btn = true
                    warning_fr = "Loading plot..."
                    spectraMz = LoadMzml(full_routeMz)
                    traceSpectra = PlotlyBase.scatter(x = mean(spectraMz[1,:]), y = mean(spectraMz[2,:]), mode="lines")
                    plotdata = [traceSpectra] # we add the data of spectra to the plot
                    spectraMz = nothing # Important for memory cleaning
                    GC.gc() # Trigger garbage collection
                    if Sys.islinux()
                        ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
                    end
                    warning_fr = "Plot loaded."
                    Disab_btn = false
                    full_routeMz2 = full_routeMz # to avoid creating the plot if its the same file read as before
                end
            else
                full_route = "/"
                warning_fr = "is not an imzML or mzML file"
            end
        catch e
            msg = "There was an error, please verify the file and try again. $(e)"
        end
    end
    @onchange Nmass begin
        indeximg = floor(Int, Nmass)
        indeximgTriq = floor(Int, Nmass)
        lastimg = floor(Int, Nmass)
        lastimgTriq = floor(Int, Nmass)
    end
    # The onbutton handler will set the variable to false after the block is executed
    @onbutton Main_Process begin
        Disab_btn = true # We disable the button to avoid multiple requests
        indeximg = floor(Int, Nmass)
        full_route = joinpath(file_route, file_name)
        if isfile(full_route) && Nmass > 0 && Tol > 0 && Tol <= 1
            msg = "File exists, Nmass=$(Nmass) Tol=$(Tol). Loading file will begin, please be patient."
            try
                spectra = LoadImzml(full_route)
                msg = "File loaded. Creating Spectra with the specific mass and tolerance, please be patient."
                slice = GetSlice(spectra, Nmass, Tol)
                fig = CairoMakie.Figure(size = (100, 200)) #container
                if triqEnabled # if we have TrIQ
                    if triqColor < 1 || triqColor > 256 ||triqProb < 0 || triqProb > 1
                        msg = "Incorrect TrIQ values, please adjust accordingly and try again."
                    else
                        SaveBitmap(joinpath("public", "TrIQ_$(floor(Int, Nmass)).bmp"),TrIQ(slice, Int(triqColor), triqProb),ViridisPalette)
                        testT = "/TrIQ_$(floor(Int, Nmass)).bmp" # we define the starting value of the images
                        msgtriq = "TrIQ image with the Nmass of $(floor(Int, Nmass))"
                        ticks = round.(range(0, stop = maximum(TrIQ(slice, Int(triqColor), triqProb)), length = 10), digits = 2)
                        Colorbar(fig[1, 1], colormap = rgb_ViridisPalette, limits = (0, maximum(TrIQ(slice, Int(triqColor), triqProb))),ticks = ticks, label = "Intensity")
                        save("public/colorbar_TrIQ_$(floor(Int, Nmass)).png", fig)
                        msg = "The file has been created successfully inside the 'public' folder of the app."
                        colorbarT = "/colorbar_TrIQ_$(floor(Int, Nmass)).png"
                    end
                else # if we don't use TrIQ
                    SaveBitmap(joinpath("public", "$(floor(Int, Nmass)).bmp"),IntQuant(slice),ViridisPalette)
                    test = "/$(floor(Int, Nmass)).bmp" # we define the starting value of the images
                    msgimg = "image with the Nmass of $(floor(Int, Nmass))"
                    ticks = round.(range(0, stop = maximum(slice), length = 10), digits = 2)
                    Colorbar(fig[1, 1], colormap = rgb_ViridisPalette, limits = (0, maximum(slice)),ticks = ticks, label = "Intensity")
                    save("public/colorbar_$(floor(Int, Nmass)).png", fig)
                    msg = "The file has been created successfully inside the 'public' folder of the app."
                    colorbar = "/colorbar_$(floor(Int, Nmass)).png"
                end
            catch e
                msg = "There was an error loading the ImzML file, please verify the file accordingly and try again. $(e)"
            end
        else
            msg = "File does not exist or a parameter is incorrect, please try again."
        end
        spectra = nothing # Important for memory cleaning
        slice = nothing
        GC.gc() # Trigger garbage collection
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
        end
        Disab_btn = false
    end
    @onbutton ImgMinus begin
        indeximg-=1
        while !isfile("public/$(indeximg).bmp") && indeximg > 0
            indeximg -= 1
        end
        if(indeximg <= 0) #if it doesn't find a lower value image
            indeximg = lastimg
        end
        test = "/$(indeximg).bmp"
        colorbar = "/colorbar_$(indeximg).png"
        msgimg = "image with the Nmass of $(indeximg)"
        lastimg = indeximg
    end
    @onbutton ImgPlus begin
        indeximg+=1
        while !isfile("public/$(indeximg).bmp") && indeximg < 2001
            indeximg += 1
        end
        if(indeximg >= 2001)
            indeximg = lastimg
        end
        test = "/$(indeximg).bmp"
        colorbar = "/colorbar_$(indeximg).png"
        msgimg = "image with the Nmass of $(indeximg)"
        lastimg = indeximg
    end

    @onbutton ImgMinusT begin
        indeximgTriq-=1
        while !isfile("public/TrIQ_$(indeximgTriq).bmp") && indeximgTriq > 0
            indeximgTriq -= 1
        end
        if(indeximgTriq <= 0) #if it doesn't find a lower value image
            indeximgTriq = lastimgTriq
        end
        testT = "/TrIQ_$(indeximgTriq).bmp"
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
        testT = "/TrIQ_$(indeximgTriq).bmp"
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
# register a new route and the page that will be loaded on access
@page("/", "app.jl.html")
end

# == Advanced features ==
#=
- The @private macro defines a reactive variable that is not sent to the browser.
This is useful for storing data that is unique to each user session but is not needed
in the UI.
    @private table = DataFrame(a = 1:10, b = 10:19, c = 20:29)

=#

