module App
# == Packages ==
# set up Genie development environment. Use the Package Manager to install new packages
#MUST FIRST USE IN TERMINAL:
#julia
#]
#add Pkg Libz
##to add the custom package go to terminal:
#add https://github.com/CINVESTAV-LABI/julia_mzML_imzML
using GenieFramework
using Pkg
using Libz
using julia_mzML_imzML 
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
    #@out test = "/wp.jpeg"#slash means it's getting the info from 'public' folder
    #@out test = "/95.bmp"
    #mkdir("new_folder")
    @in file_route = ""
    @in file_name = ""
    @out warning_fr = ""
    @in Nmass = 0.0
    @in Tol = 0.0
    @in triqEnabled = false
    @in triqProb = 0.0
    @in Main_Process = false
    @out indeximg = 0
    @out lastimg = 0
    @in ImgPlus = false
    @in ImgMinus = false
    @in msg = ""
    @out full_route = ""

    # == Reactive handlers ==
    # reactive handlers watch a variable and execute a block of code when
    # its value changes
    @onchange triqEnabled begin
        if !triqEnabled
            triqProb = 0.0
        end
    end
    @onchange file_name begin
        if contains_word(file_name, ".imzML")
            warning_fr = ""
            full_route = joinpath( file_route, file_name )
        else
            warning_fr = "is not an imzML file"
        end
    end
    @onchange Nmass begin
        indeximg = Nmass
        lastimg = Nmass
    end
    
    # the onbutton handler will set the variable to false after the block is executed
    #/home/julian/Documentos/Cinvestav_2024/Web/Archivos IMZML/royaimg.imzML
    @onbutton Main_Process begin
        indeximg = Nmass
        full_route = joinpath( file_route, file_name )
        if isfile(full_route)
            msg="File exists, Nmass=$(Nmass) Tol=$(Tol). Please do not press the start button until confirmation"
            spectra = LoadImzml(full_route)
            msg = "File loaded. Please do not press the start button until confirmation"
            slice = GetSlice(spectra, Nmass, Tol)
            if triqProb != 0    
                # SaveBitmap( joinpath( file_route, "TrIQ_$(Int(Nmass)).bmp" ),
                # TrIQ( slice, Int(Nmass), triqProb ),
                # ViridisPalette )
                SaveBitmap( joinpath( "public", "TrIQ_$(Int(Nmass)).bmp" ),
                TrIQ( slice, Int(Nmass), triqProb ),
                ViridisPalette )
            else
                # SaveBitmap( joinpath( file_route, "$(Int(Nmass)).bmp" ),
                # IntQuant( slice ),
                # ViridisPalette )
                #test
                SaveBitmap( joinpath( "public", "$(Int(Nmass)).bmp" ),
                IntQuant( slice ),
                ViridisPalette )
            end
            msg="The file has been created inside the 'public' folder of the app"
        else
            msg="File does not exist"
        end
        spectra = nothing   
        slice = nothing
        test = "/$(Int(Nmass)).bmp"#we define the starting value of the images
        GC.gc()
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
        msg = "image with the Nmass of $(indeximg)"
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
        msg = "image with the Nmass of $(indeximg)"
        lastimg = indeximg
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