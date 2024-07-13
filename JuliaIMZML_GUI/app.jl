module App
# == Packages ==
# set up Genie development environment. Use the Package Manager to install new packages
#MUST FIRST USE IN TERMINAL:
#julia
#]
#add Pkg Libz
using GenieFramework
using Pkg
using Libz
#to add the custom package go to terminal:
#julia
#]
#add https://github.com/CINVESTAV-LABI/julia_mzML_imzML
#Pkg.activate("/julia_mzML_imzML")
using julia_mzML_imzML #could  add this to a try-catch to verify they're putting the correct name of the file
@genietools

# == Code import ==
# add your data analysis code here or in the lib folder. Code in lib/ will be
# automatically loaded 
function mean_value(x)
    sum(x) / length(x)
end
function contains_word(x, word)
    occursin(word, x)
end

# == Reactive code ==
# add reactive code to make the UI interactive
@app begin
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block
    @in N = 10
    @out x = collect(1:10)
    @out y = randn(10) # plot data must be an array or a DataFrame
    @out msg = "The average is 0."
    @in shift = false

    #@out test = "/Borrador_Julia_GUI.jpeg"
    @out test = "/wp"
    #mkdir("new_folder")
    @in file_route = "None"
    @out warning_fr = ""
    @in Nmass = 0.0
    @in Tol = 0.0
    @in triqEnabled = false
    @in triqProb = 0.0
    @in Main_Process = false
    @out indeximg = 0
    @in ImgPlus = false
    @in ImgMinus = false

    @onchange triqEnabled begin
        if !triqEnabled
            triqProb = 0.0
        end
    end
    @onchange file_route begin
        if contains_word(file_route, ".imzML")
            warning_fr = ""
        else
            warning_fr = "is not an imzML file"
        end
    end
    #/home/julian/Documentos/Cinvestav_2024/Web/Archivos IMZML/royaimg.imzML
    @onbutton Main_Process begin
        test = "/Borrador_Julia_GUI.jpeg"
        indeximg = Nmass
        #spectra = LoadImzml(file_route)
        #slice = getSlice(spectra, Nmass, Tol)
        #SaveBitmap( "$(Nmass).bmp",
        #IntQuant( slice ),
        #ViridisPalette )
        samplesDir = "/home/julian/Documentos/Cinvestav_2024/Web/Archivos IMZML"
        fileName = "royaimg.imzML"
        if isfile("/home/julian/Documentos/Cinvestav_2024/Web/Archivos IMZML/royaimg.imzML")
            msg="File exists, Nmass=$(Nmass) Tol=$(Tol)"
        else
            msg="File does not exist"
        end
        spectra = LoadImzml(joinpath( samplesDir, fileName ))
        slice = getSlice(spectra, Nmass, Tol)
        SaveBitmap( joinpath( samplesDir, "$(Nmass).bmp" ),
        IntQuant( slice ),
        ViridisPalette )
    end
    @onbutton ImgMinus begin
        indeximg-=1
        while !isfile("/$(indeximg).jpeg") && indeximg > 0 && indeximg <999
            indeximg-=1
        end
        test = "/$(indeximg).jpeg"
        msg = string(indeximg)
    end
    @onbutton ImgPlus begin
        indeximg+=1
        while !isfile("/$(indeximg).jpeg") && indeximg > 0 && indeximg <999
            indeximg+=1
        end
        test = "/$(indeximg).jpeg"
        msg = string(indeximg)
    end

    # == Reactive handlers ==
    # reactive handlers watch a variable and execute a block of code when
    # its value changes
    @onchange N begin
        # the values of x, result and msg in the UI will
        # be automatically updated
        x = collect(1:N)
        y = rand(N)
        result = mean_value(rand(N))
        msg = "The average is $result."
    end
    # the onbutton handler will set the variable to false after the block is executed
    @onbutton shift begin
        y = circshift(y, 1)
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