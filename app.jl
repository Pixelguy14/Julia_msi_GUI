# app.jl

module App
# ==Packages ==
using GenieFramework # Set up Genie development environment.
using Pkg
using Libz
using PlotlyBase
using CairoMakie
using Colors
using MSI_src # Import the new MSIData library
using Statistics
using NaturalSort
using Images
using LinearAlgebra
using NativeFileDialog # Opens the file explorer depending on the OS
using StipplePlotly
using Base.Filesystem: mv # To rename files in the system
using Printf # Required for @sprintf macro in colorbar generation
using JSON
using Dates
using Base.Threads

# Bring MSIData into App module's scope
using .MSI_src: MSIData, OpenMSIData, process_spectrum, IterateSpectra, ImzMLSource, _iterate_spectra_fast, MzMLSource, find_mass, ViridisPalette, get_mz_slice, get_multiple_mz_slices, quantize_intensity, save_bitmap, median_filter, save_bitmap, downsample_spectrum, TrIQ, precompute_analytics, ImportMzmlFile, generate_colorbar_image, load_and_prepare_mask, set_global_mz_range!

if !@isdefined(increment_image)
    include("./julia_imzML_visual.jl")
end

# --- Memory Validation Logging ---
if get(ENV, "GENIE_ENV", "dev") != "prod"
    function get_rss_mb()
        if !Sys.islinux()
            return 0.0
        end
        try
            pid = getpid()
            cmd = `ps -p $pid -o rss=`
            rss_kb_str = read(cmd, String)
            rss_kb = parse(Int, strip(rss_kb_str))
            return round(rss_kb / 1024, digits=2)
        catch e
            @warn "Could not get RSS via `ps` command. Error: $e"
            return 0.0
        end
    end

    function log_memory_usage(context::String, msi_data_val)
        rss_mb = get_rss_mb()
        
        msi_data_size_mb = 0
        if msi_data_val !== nothing
            msi_data_size_mb = round(Base.summarysize(msi_data_val) / (1024^2), digits=2)
        end

        gc_time_s = round(GC.time(), digits=3)
        
        println("--- MEMORY LOG [$(context)] ---")
        println("  Timestamp: $(now())")
        println("  Process RSS: $(rss_mb) MB")
        println("  msi_data size: $(msi_data_size_mb) MB")
        println("  Cumulative GC time: $(gc_time_s) s")
        println("--------------------------")
    end
else
    log_memory_usage(context::String, msi_data_val) = nothing # No-op for production
end

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
    @out btnSpectraDisable=false
    # Loading animations
    @in progress=false
    @in progressPlot=false
    @in progressSpectraPlot=false
    # Text field validations
    @in triqEnabled=false
    @in SpectraEnabled=false
    @in MFilterEnabled=false
    @in maskEnabled=false
    # Dialogs
    @in warning_msg=false
    @in CompareDialog=false

    ## Interface Variables
    @in file_route=""
    @in file_name=""
    @in Nmass="0.0"
    @in Tol=0.1
    @in triqProb=0.98
    @in colorLevel=20

    ## Interface Buttons
    @in btnSearch=false # To search for files in your device
    @in btnAddBatch = false
    @in clear_batch_btn = false
    @out batch_file_count = 0
    @in mainProcess=false # To generate images
    @in compareBtn=false # To open dialog
    @in createMeanPlot=false # To generate mean spectrum plot
    @in createXYPlot=false # To generate an spectrum plot according to the xy values inputed
    @in createSumPlot=false # To generate a sum of all the spectrum plots
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
    @in imgPlusCompLeft=false
    @in imgMinusCompLeft=false
    @in imgPlusTCompLeft=false
    @in imgMinusTCompLeft=false
    @in imgPlusCompRight=false
    @in imgMinusCompRight=false
    @in imgPlusTCompRight=false
    @in imgMinusTCompRight=false

    ## Tabulation variables
    @out tabIDs=["tab0","tab1","tab2","tab3","tab4"]
    @out tabLabels=["Image", "TrIQ", "Spectrum Plot", "Topography Plot","Surface Plot"]
    @in selectedTab="tab0"
    @out CompTabIDsLeft=["tab0","tab1","tab2","tab3","tab4"]
    @out CompTabLabelsLeft=["Image", "TrIQ", "Spectrum Plot", "Topography Plot","Surface Plot"]
    @in CompSelectedTabLeft="tab0"
    @out CompTabIDsRight=["tab0","tab1","tab2","tab3","tab4"]
    @out CompTabLabelsRight=["Image", "TrIQ", "Spectrum Plot", "Topography Plot","Surface Plot"]
    @in CompSelectedTabRight="tab0"

    # Interface Images
    @out imgInt="/.bmp" # image Interface
    @out imgIntT="/.bmp" # image Interface TrIQ
    @out colorbar="/.png"
    @out colorbarT="/.png"
    # Interface controlling for the comparative view
    @out imgIntCompLeft="/.bmp"
    @out imgIntTCompLeft="/.bmp"
    @out colorbarCompLeft="/.png"
    @out colorbarTCompLeft="/.png"
    @out imgIntCompRight="/.bmp"
    @out imgIntTCompRight="/.bmp"
    @out colorbarCompRight="/.png"
    @out colorbarTCompRight="/.png"
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
    @out msgimgCompLeft=""
    @out msgtriqCompLeft=""
    @out msgimgCompRight=""
    @out msgtriqCompRight=""

    # Centralized MSIData object
    @out msi_data::Union{MSIData, Nothing} = nothing

    # Metadata table variables
    @in showMetadataDialog = false
    @in showMetadataBtn = false
    @out metadata_columns = []
    @out metadata_rows = []
    @out btnMetadataDisable = false
    @in selected_folder_metadata = ""

    # Saves the route where imzML and mzML files are located
    @out full_route=""

    # == Converter Tab Variables ==
    @in left_tab = "generator"
    @out mzml_full_route = ""
    @out sync_full_route = ""
    @in btnSearchMzml = false
    @in btnSearchSync = false
    @in convert_process = false
    @out progress_conversion = false
    @out msg_conversion = ""
    @out btnConvertDisable = true

    # == Pre Processing Variables ==
    @in pre_tab = "stabilization"

    # == Batch Summary Dialog ==
    @in showBatchSummary = false
    @out batch_summary = ""

    # == Batch Processing & Registry Variables ==
    @private registry_init_done = false
    @in selected_files = String[]
    @in available_folders = String[]
    @in image_available_folders = String[]
    @out registry_path = abspath(joinpath(@__DIR__, "public", "registry.json"))
    # Progress reporting
    @out overall_progress = 0.0
    @out progress_message = ""

    # == Folder-based UI State ==
    @in selected_folder_main = ""
    @in selected_folder_compare_left = ""
    @in selected_folder_compare_right = ""


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
    @out current_msiCompLeft=""
    @out current_col_msiCompLeft=""
    @out current_triqCompLeft=""
    @out current_col_triqCompLeft=""
    @out current_msiCompRight=""
    @out current_col_msiCompRight=""
    @out current_triqCompRight=""
    @out current_col_triqCompRight=""

    ## Time measurement variables
    @out sTime=time()
    @out fTime=time()
    @out eTime=time()

    ## Plots
    # Local image to plot
    layoutImg=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=14,
                color="black"
            )
        ),
        xaxis=PlotlyBase.attr(
            visible=false,
            scaleanchor="y",
            range=[0, 0]
        ),
        yaxis=PlotlyBase.attr(
            visible=false,
            range=[0, 0]
        ),
        margin=attr(l=0,r=0,t=0,b=0,pad=0)
    )
    traceImg=PlotlyBase.heatmap(x=Vector{Float64}(), y=Vector{Float64}())
    @out plotdataImg=[traceImg]
    @out plotlayoutImg=layoutImg
    # For the image in the comparative view
    @out plotdataImgCompLeft=[traceImg]
    @out plotlayoutImgCompLeft=layoutImg
    @out plotdataImgCompRight=[traceImg]
    @out plotlayoutImgCompRight=layoutImg

    # For triq image 
    @out plotdataImgT=[traceImg]
    @out plotlayoutImgT=layoutImg
    # For the triq image in the comparative view
    @out plotdataImgTCompLeft=[traceImg]
    @out plotlayoutImgTCompLeft=layoutImg
    @out plotdataImgTCompRight=[traceImg]
    @out plotlayoutImgTCompRight=layoutImg
    # Interface Plot Spectrum
    layoutSpectra=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="Spectrum plot",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        hovermode="closest",
        xaxis=PlotlyBase.attr(
            title="<i>m/z</i>",
            showgrid=true
        ),
        yaxis=PlotlyBase.attr(
            title="Intensity",
            showgrid=true,
            tickformat = ".3g"
        ),
        margin=attr(l=0,r=0,t=120,b=0,pad=0)
    )
    # Dummy 2D scatter plot
    traceSpectra=PlotlyBase.scatter(x=Vector{Float64}(), y=Vector{Float64}(), mode="lines", marker=attr(size=1, color="blue", opacity=0.1))
    # Create conection to frontend
    @out plotdata=[traceSpectra]
    @out plotlayout=layoutSpectra
    @in xCoord=0
    @in yCoord=0
    @out xSpectraMz = Vector{Float64}()
    @out ySpectraMz = Vector{Float64}()

    # Interactive plot reactions
    @in data_click=Dict{String,Any}()
    #@in data_selected=Dict{String,Any}() # Selected is for areas, this can work for the masks
    #<plotly id="plotStyle" :data="plotdata" :layout="plotlayout" @click="data_selected" class="q-pa-none q-ma-none sync_data"></plotly>

    # Interface Plot Surface
    layoutContour=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="2D Topographic map",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
        xaxis=PlotlyBase.attr(
            visible=false,
            scaleanchor="y"
        ),
        yaxis=PlotlyBase.attr(
            visible=false
        ),
        margin=attr(l=0,r=0,t=100,b=0,pad=0)
    )
    # Dummy 2D surface plot
    traceContour=PlotlyBase.contour(x=Vector{Float64}(), y=Vector{Float64}(), mode="lines")
    # Create conection to frontend
    @out plotdataC=[traceContour]
    @out plotlayoutC=layoutContour

    # Interface Plot 3d
    # Define the layout for the 3D plot
    layout3D=PlotlyBase.Layout(
        title=PlotlyBase.attr(
            text="3D Surface plot",
            font=PlotlyBase.attr(
                family="Roboto, Lato, sans-serif",
                size=18,
                color="black"
            )
        ),
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
    trace3D=PlotlyBase.surface(x=Vector{Float64}(), y=Vector{Float64}(), z=Matrix{Float64}(undef, 0, 0),
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

    # This handler correctly uses pick_file and loads the selected file
    # as the active dataset for the UI.
    @onbutton btnSearch begin
        picked_route = pick_file(; filterlist="imzML,imzml,mzML,mzml")
        if isempty(picked_route)
            return
        end

        progress = true
        msg = "Opening file: $(basename(picked_route))..."

        try
            dataset_name = replace(basename(picked_route), r"(\.(imzML|imzml|mzML|mzml))$"i => "")
            registry = load_registry(registry_path)
            existing_entry = get(registry, dataset_name, nothing)

            # --- Fast Load Path ---
            is_same_file = (existing_entry !== nothing && existing_entry["source_path"] == picked_route)
            if is_same_file && !isempty(get(existing_entry, "metadata", Dict()))
                msg = "Fast loading pre-processed file: $(dataset_name)"
                println(msg)
                
                full_route = existing_entry["source_path"]
                metadata_rows = existing_entry["metadata"]["summary"]
                
                dims_str = first(filter(r -> r["parameter"] == "Image Dimensions", metadata_rows))["value"]
                dims = parse.(Int, split(dims_str, " x "))
                imgWidth, imgHeight = dims[1], dims[2]

                msi_data = nothing # Ensure data is not held in memory
                log_memory_usage("Fast Load (msi_data cleared)", msi_data)
                btnMetadataDisable = false
                btnStartDisable = false
                btnPlotDisable = false
                btnSpectraDisable = false
                SpectraEnabled = true
                selected_folder_main = dataset_name
                
                # Update folder lists in UI
                all_folders = sort(collect(keys(registry)), lt=natural)
                img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)
                available_folders = deepcopy(all_folders)
                image_available_folders = deepcopy(img_folders)

                msg = "Successfully loaded pre-processed dataset: $(dataset_name)"
                progress = false
                return
            end

            # --- Full Load Path ---
            msg = "Performing first-time analysis for: $(basename(picked_route))..."
            local local_full_route
            if endswith(picked_route, r"imzml"i)
                local_full_route = replace(picked_route, r"\.imzml$"i => ".imzML")
                if picked_route != local_full_route
                    mv(picked_route, local_full_route, force=true)
                end
            else
                local_full_route = picked_route
            end
            full_route = local_full_route

            sTime = time()
            loaded_data = OpenMSIData(local_full_route)
            is_imzML = loaded_data.source isa ImzMLSource
            
            precompute_analytics(loaded_data)
            
            metadata_columns = [
                Dict("name" => "parameter", "label" => "Parameter", "field" => "parameter", "align" => "left"),
                Dict("name" => "value", "label" => "Value", "field" => "value", "align" => "left"),
            ]
            summary_stats = extract_metadata(loaded_data, local_full_route)
            metadata_rows = summary_stats["summary"]
            btnMetadataDisable = isempty(metadata_rows)

            w, h = loaded_data.image_dims
            imgWidth, imgHeight = w > 0 ? (w, h) : (500, 500)

            update_registry(registry_path, dataset_name, local_full_route, summary_stats, is_imzML)
            
            # Update folder lists in UI
            registry = load_registry(registry_path)
            all_folders = sort(collect(keys(registry)), lt=natural)
            img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)
            available_folders = deepcopy(all_folders)
            image_available_folders = deepcopy(img_folders)

            selected_folder_main = dataset_name
            msi_data = loaded_data
            log_memory_usage("Full Load", msi_data)

            eTime = round(time() - sTime, digits=3)
            msg = "Active file loaded in $(eTime) seconds. Dataset '$(dataset_name)' is ready for analysis."

            btnStartDisable = false
            btnPlotDisable = false
            btnSpectraDisable = false
            SpectraEnabled = true
            
        catch e
            msi_data = nothing
            msg = "Error loading active file: $e"
            warning_msg = true
            btnStartDisable = true
            btnSpectraDisable = true
            SpectraEnabled = false
            btnMetadataDisable = true
            @error "File loading failed" exception=(e, catch_backtrace())
        finally
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
            progress = false
            progressSpectraPlot = false
        end
    end

    # This new handler correctly adds the file from full_route to the batch list.
    @onbutton btnAddBatch begin
        if isempty(full_route) || full_route == "unknown (manually added)"
            msg = "No active file selected to add to batch."
            warning_msg = true
            return
        end

        if !(full_route in selected_files)
            push!(selected_files, full_route)
            selected_files = deepcopy(selected_files) # Force reactivity
            batch_file_count = length(selected_files)
            msg = "File added to batch."
        else
            msg = "File is already in the batch list."
            warning_msg = true
        end
    end

    @onbutton clear_batch_btn begin
        selected_files = String[]
        batch_file_count = 0
        msg = "Batch cleared"
    end

    @onchange selected_files begin
        batch_file_count = length(selected_files)
    end

    @onchange full_route begin
        if !isempty(full_route) && !(full_route in selected_files)
            push!(selected_files, full_route)
            selected_files = deepcopy(selected_files) # Force reactivity
            batch_file_count = length(selected_files)
            msg = "File automatically added to batch"
        end
    end

    @onbutton showMetadataBtn begin
        if !isempty(available_folders)
            if !isempty(selected_folder_main)
                selected_folder_metadata = selected_folder_main
            elseif !isempty(available_folders)
                selected_folder_metadata = first(available_folders)
            end
            showMetadataDialog = true
        else
            msg = "No processed datasets available."
            warning_msg = true
        end
    end

    @onchange selected_folder_metadata begin
        if !isempty(selected_folder_metadata)
            registry = load_registry(registry_path)
            dataset_info = get(registry, selected_folder_metadata, nothing)

            if dataset_info !== nothing && haskey(dataset_info, "metadata") && !isempty(get(dataset_info["metadata"], "summary", []))
                metadata_rows = dataset_info["metadata"]["summary"]
                btnMetadataDisable = false
            else
                metadata_rows = []
                btnMetadataDisable = true
                msg = "Metadata not found in registry for $(selected_folder_metadata)."
            end
        end
    end

    @onchange btnSearchMzml, btnSearchSync begin
        if btnSearchMzml
            picked_route = pick_file(; filterlist="mzML,mzml")
            if !isempty(picked_route)
                mzml_full_route = picked_route
            end
            btnSearchMzml = false # Reset the button
        end

        if btnSearchSync
            picked_route = pick_file(; filterlist="txt")
            if !isempty(picked_route)
                sync_full_route = picked_route
            end
            btnSearchSync = false # Reset the button
        end

        # Enable button only if both files are selected
        btnConvertDisable = isempty(mzml_full_route) || isempty(sync_full_route)
    end

    @onbutton convert_process begin
        if isempty(mzml_full_route) || isempty(sync_full_route)
            msg_conversion = "Please select both an .mzML file and a .txt sync file."
            warning_msg = true
            return
        end

        progress_conversion = true
        btnConvertDisable = true
        msg_conversion = "Starting conversion process..."

        try
            sTime = time()
            target_imzml = replace(mzml_full_route, r"\.(mzml|mzML)$" => ".imzML")
            
            msg_conversion = "Converting $(basename(mzml_full_route)) to $(basename(target_imzml))... This may take a while."

            success = ImportMzmlFile(mzml_full_route, sync_full_route, target_imzml)

            fTime = time()
            eTime = round(fTime - sTime, digits=3)

            if success
                msg_conversion = "Conversion successful in $(eTime) seconds. Output file: $(basename(target_imzml))"
            else
                msg_conversion = "Conversion failed after $(eTime) seconds. Check console for errors."
                warning_msg = true
            end

        catch e
            msg_conversion = "An error occurred during conversion: $e"
            warning_msg = true
            @error "Conversion failed" exception=(e, catch_backtrace())
        finally
            progress_conversion = false
            # Re-enable button if files are still selected
            btnConvertDisable = isempty(mzml_full_route) || isempty(sync_full_route)
        end
    end
    
    @onbutton mainProcess @time begin
        # --- UI State Update ---
        progress = true
        btnStartDisable = true
        btnPlotDisable = true
        btnSpectraDisable = true
        overall_progress = 0.0
        progress_message = "Preparing batch process..."

        # --- CAPTURE CURRENT VALUES HERE ---
        current_selected_files = selected_files
        current_nmass = Nmass
        current_tol = Tol
        current_color_level = colorLevel
        current_triq_enabled = triqEnabled
        current_triq_prob = triqProb
        current_mfilter_enabled = MFilterEnabled
        current_mask_enabled = maskEnabled
        current_registry_path = registry_path

        println("starting main process with $(length(current_selected_files)) files")
        total_time_start = time()
        try
            # --- 1. Parameter Validation ---
            if isempty(current_selected_files)
                progress_message = "No .imzML files in batch. Please add files first."
                warning_msg = true
                println(progress_message)
                return
            end

            masses = Float64[]
            try
                masses = [parse(Float64, strip(m)) for m in split(current_nmass, ',', keepempty=false)]
            catch e
                progress_message = "Invalid m/z value(s). Please provide a comma-separated list of numbers. Error: $e"
                warning_msg = true
                return
            end

            if isempty(masses)
                progress_message = "No valid m/z values found. Please provide comma-separated positive numbers."
                warning_msg = true
                return
            end

            # --- 2. Batch Processing Loop ---
            num_files = length(current_selected_files)
            total_steps = num_files
            current_step = 0
            errors = Dict("load_errors" => String[], "slice_errors" => String[], "io_errors" => String[])
            newly_created_folders = String[]
            files_without_mask = 0

            for (file_idx, file_path) in enumerate(current_selected_files)
                progress_message = "Processing file $file_idx/$num_files: $(basename(file_path))"
                overall_progress = current_step / total_steps
                
                all_params = (
                    tolerance = current_tol,
                    colorL = current_color_level,
                    triqE = current_triq_enabled,
                    triqP = current_triq_prob,
                    medianF = current_mfilter_enabled,
                    registry = current_registry_path,
                    fileIdx = file_idx,
                    nFiles = num_files
                )

                success, error_msg = process_file_safely(file_path, masses, all_params, progress_message, overall_progress, use_mask=current_mask_enabled)

                if !success
                    push!(errors["load_errors"], error_msg)
                else
                    push!(newly_created_folders, replace(basename(file_path), r"\.imzML$"i => ""))
                end
                current_step += 1
            end

            # --- 3. Final Report ---
            total_time_end = round(time() - total_time_start, digits=3)
            
            registry = load_registry(current_registry_path)
            all_folders = sort(collect(keys(registry)), lt=natural)
            img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)
            available_folders = deepcopy(all_folders)
            image_available_folders = deepcopy(img_folders)

            if !isempty(newly_created_folders)
                selected_folder_main = first(newly_created_folders)
            end

            successful_files = length(newly_created_folders)
            total_errors = sum(length, values(errors))

            if total_errors == 0
                msg = "Successfully processed all $(successful_files) file(s) in $(total_time_end) seconds."
            else
                msg = "Batch completed in $(total_time_end) seconds with $(total_errors) error(s)."
                warning_msg = true
            end

            mask_summary = current_mask_enabled ? "\nFiles processed without a mask: $(files_without_mask)" : ""

            batch_summary = """
            Processed $(successful_files)/$(num_files) files successfully.
            $(mask_summary)

            Errors by category:
            • Load failures: $(length(errors["load_errors"]))
            • Slice generation: $(length(errors["slice_errors"]))  
            • I/O issues: $(length(errors["io_errors"]))

            Detailed errors:
            $(join(vcat(values(errors)...), "\n"))
            """
            showBatchSummary = true

            # Update UI to display the last generated image
            if !isempty(newly_created_folders)
                timestamp = string(time_ns())
                folder_path = joinpath("public", selected_folder_main)

                if current_triq_enabled
                    triq_files = filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path))
                    col_triq_files = filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path))

                    if !isempty(triq_files)
                        latest_triq = triq_files[argmax([mtime(joinpath(folder_path, f)) for f in triq_files])]
                        current_triq = latest_triq
                        imgIntT = "/$(selected_folder_main)/$(current_triq)?t=$(timestamp)"
                        plotdataImgT, plotlayoutImgT, _, _ = loadImgPlot(imgIntT)
                        text_nmass = replace(current_triq, r"TrIQ_|.bmp" => "")
                        msgtriq = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
                        
                        if !isempty(col_triq_files)
                            latest_col_triq = col_triq_files[argmax([mtime(joinpath(folder_path, f)) for f in col_triq_files])]
                            current_col_triq = latest_col_triq
                            colorbarT = "/$(selected_folder_main)/$(current_col_triq)?t=$(timestamp)"
                        else
                            colorbarT = ""
                        end
                        selectedTab = "tab1"
                    end
                else # Not TrIQ enabled, display regular MSI image
                    msi_files = filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path))
                    col_msi_files = filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path))

                    if !isempty(msi_files)
                        latest_msi = msi_files[argmax([mtime(joinpath(folder_path, f)) for f in msi_files])]
                        current_msi = latest_msi
                        imgInt = "/$(selected_folder_main)/$(current_msi)?t=$(timestamp)"
                        plotdataImg, plotlayoutImg, _, _ = loadImgPlot(imgInt)
                        text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
                        msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"

                        if !isempty(col_msi_files)
                            latest_col_msi = col_msi_files[argmax([mtime(joinpath(folder_path, f)) for f in col_msi_files])]
                            current_col_msi = latest_col_msi
                            colorbar = "/$(selected_folder_main)/$(current_col_msi)?t=$(timestamp)"
                        else
                            colorbar = ""
                        end
                        selectedTab = "tab0"
                    end
                end
            end

        catch e
            println("Error in main process: $e")
            msg = "Batch processing failed: $e"
            warning_msg = true
            @error "Main process failed" exception=(e, catch_backtrace())
        finally
            # --- UI State Reset ---
            progress = false
            btnStartDisable = false
            btnPlotDisable = false
            btnOpticalDisable = false
            btnSpectraDisable = false
            SpectraEnabled = true
            overall_progress = 0.0
            println("Done")
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end

    @onbutton createMeanPlot @time begin
        if isempty(selected_folder_main)
            msg = "No dataset selected. Please process a file and select a folder first."
            warning_msg = true
            return
        end

        progressSpectraPlot = true
        btnPlotDisable = true
        btnStartDisable = true
        msg = "Loading plot for $(selected_folder_main)..."

        try
            sTime = time()
            registry = load_registry(registry_path)
            entry = registry[selected_folder_main]
            target_path = entry["source_path"]

            if target_path == "unknown (manually added)"
                msg = "Dataset selected contained no route."
                warning_msg = true
                return
            end

            if msi_data === nothing || full_route != target_path
                if msi_data !== nothing
                    close(msi_data)
                end
                msg = "Reloading $(basename(target_path)) for analysis..."
                full_route = target_path
                msi_data = OpenMSIData(target_path)
                if haskey(get(entry, "metadata", Dict()), "global_min_mz") && entry["metadata"]["global_min_mz"] !== nothing
                    raw_min = entry["metadata"]["global_min_mz"]
                    raw_max = entry["metadata"]["global_max_mz"]
                    min_val = isa(raw_min, Dict) ? get(raw_min, "value", raw_min) : raw_min
                    max_val = isa(raw_max, Dict) ? get(raw_max, "value", raw_max) : raw_max
                    set_global_mz_range!(msi_data, convert(Float64, min_val), convert(Float64, max_val))
                else
                    precompute_analytics(msi_data)
                end
            end

            local mask_path_for_plot::Union{String, Nothing} = nothing
            if maskEnabled && get(entry, "has_mask", false)
                mask_path_for_plot = get(entry, "mask_path", "")
                if !isfile(mask_path_for_plot)
                    @warn "Mask not found for plotting: $(mask_path_for_plot). Plotting without mask."
                    mask_path_for_plot = nothing
                end
            end

            plotdata, plotlayout, xSpectraMz, ySpectraMz = meanSpectrumPlot(msi_data, selected_folder_main, mask_path=mask_path_for_plot)
            selectedTab = "tab2"
            fTime = time()
            eTime = round(fTime - sTime, digits=3)
            msg = "Plot loaded in $(eTime) seconds"
            log_memory_usage("Mean Plot Generated", msi_data)
        catch e
            msg = "Could not generate mean spectrum plot: $e"
            warning_msg = true
            @error "Mean spectrum plotting failed" exception=(e, catch_backtrace())
        finally
            progressSpectraPlot = false
            btnPlotDisable = false
            btnSpectraDisable = false
            btnStartDisable = false
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end

    @onbutton createSumPlot @time begin
        if isempty(selected_folder_main)
            msg = "No dataset selected. Please process a file and select a folder first."
            warning_msg = true
            return
        end

        progressSpectraPlot = true
        btnPlotDisable = true
        btnStartDisable = true
                    msg = "Loading total spectrum plot for $(selected_folder_main)..."
        
                    try
                        sTime = time()
                        registry = load_registry(registry_path)
                        entry = registry[selected_folder_main]
                        target_path = entry["source_path"]
        
                        if target_path == "unknown (manually added)"
                            msg = "Dataset selected contained no route."
                            warning_msg = true
                            return
                        end
        
                        if msi_data === nothing || full_route != target_path
                            if msi_data !== nothing
                                close(msi_data)
                            end
                            msg = "Reloading $(basename(target_path)) for analysis..."
                            full_route = target_path
                            msi_data = OpenMSIData(target_path)
                            if haskey(get(entry, "metadata", Dict()), "global_min_mz") && entry["metadata"]["global_min_mz"] !== nothing
                                raw_min = entry["metadata"]["global_min_mz"]
                                raw_max = entry["metadata"]["global_max_mz"]
                                min_val = isa(raw_min, Dict) ? get(raw_min, "value", raw_min) : raw_min
                                max_val = isa(raw_max, Dict) ? get(raw_max, "value", raw_max) : raw_max
                                set_global_mz_range!(msi_data, convert(Float64, min_val), convert(Float64, max_val))
                            else
                                precompute_analytics(msi_data)
                            end
                        end
            local mask_path_for_plot::Union{String, Nothing} = nothing
            if maskEnabled && get(entry, "has_mask", false)
                mask_path_for_plot = get(entry, "mask_path", "")
                if !isfile(mask_path_for_plot)
                    @warn "Mask not found for plotting: $(mask_path_for_plot). Plotting without mask."
                    mask_path_for_plot = nothing
                end
            end

            plotdata, plotlayout, xSpectraMz, ySpectraMz = sumSpectrumPlot(msi_data, selected_folder_main, mask_path=mask_path_for_plot)
            selectedTab = "tab2"
            fTime = time()
            eTime = round(fTime - sTime, digits=3)
            msg = "Total plot loaded in $(eTime) seconds"
            log_memory_usage("Sum Plot Generated", msi_data)
        catch e
            msg = "Could not generate total spectrum plot: $e"
            warning_msg = true
            @error "Total spectrum plotting failed" exception=(e, catch_backtrace())
        finally
            progressSpectraPlot = false
            btnPlotDisable = false
            btnSpectraDisable = false
            btnStartDisable = false
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end

    @onbutton createXYPlot @time begin
        if isempty(selected_folder_main)
            msg = "No dataset selected. Please process a file and select a folder first."
            warning_msg = true
            return
        end

        progressSpectraPlot = true
        btnStartDisable = true
        btnPlotDisable = true
        btnSpectraDisable = true
        msg = "Loading plot for $(selected_folder_main)..."

        try
            sTime = time()
            registry = load_registry(registry_path)
            
            # Add error handling for registry access
            if !haskey(registry, selected_folder_main)
                msg = "Dataset '$selected_folder_main' not found in registry."
                warning_msg = true
                return
            end
            
            entry = registry[selected_folder_main]
            target_path = entry["source_path"]

            if target_path == "unknown (manually added)"
                msg = "Dataset selected contained no route."
                warning_msg = true
                return
            end

            if msi_data === nothing || full_route != target_path
                if msi_data !== nothing
                    close(msi_data)
                end
                msg = "Reloading $(basename(target_path)) for analysis..."
                full_route = target_path
                msi_data = OpenMSIData(target_path)
                if haskey(get(entry, "metadata", Dict()), "global_min_mz") && entry["metadata"]["global_min_mz"] !== nothing
                    raw_min = entry["metadata"]["global_min_mz"]
                    raw_max = entry["metadata"]["global_max_mz"]
                    min_val = isa(raw_min, Dict) ? get(raw_min, "value", raw_min) : raw_min
                    max_val = isa(raw_max, Dict) ? get(raw_max, "value", raw_max) : raw_max
                    set_global_mz_range!(msi_data, convert(Float64, min_val), convert(Float64, max_val))
                else
                    precompute_analytics(msi_data)
                end
            end

            local mask_path_for_plot::Union{String, Nothing} = nothing
            if maskEnabled && get(entry, "has_mask", false)
                mask_path_for_plot = get(entry, "mask_path", "")
                if !isfile(mask_path_for_plot)
                    @warn "Mask not found for plotting: $(mask_path_for_plot). Plotting without mask."
                    mask_path_for_plot = nothing
                end
            end

            # Convert to positive coordinates for processing
            y_positive = yCoord < 0 ? abs(yCoord) : yCoord
            plotdata, plotlayout, xSpectraMz, ySpectraMz = xySpectrumPlot(msi_data, xCoord, y_positive, imgWidth, imgHeight, selected_folder_main, mask_path=mask_path_for_plot)
            
            # Update coordinates based on actual plot title
            # Extract title text from the Dict safely
            actual_title = if plotlayout.title isa Dict && haskey(plotlayout.title, :text)
                plotlayout.title[:text]
            elseif plotlayout.title isa Dict && haskey(plotlayout.title, "text")
                plotlayout.title["text"]
            else
                string(plotlayout.title)  # Fallback
            end
            
            if occursin("Masked Spectrum at", actual_title)
                # Extract coordinates from masked spectrum title
                coords_match = match(r"Masked Spectrum at \((\d+), (\d+)\)", actual_title)
                if coords_match !== nothing
                    xCoord = parse(Int, coords_match.captures[1])
                    yCoord = -parse(Int, coords_match.captures[2])  # Negative for display
                end
            elseif occursin("Spectrum at", actual_title)
                # Extract coordinates from regular spectrum title
                coords_match = match(r"Spectrum at \((\d+), (\d+)\)", actual_title)
                if coords_match !== nothing
                    xCoord = parse(Int, coords_match.captures[1])
                    yCoord = -parse(Int, coords_match.captures[2])  # Negative for display
                end
            else
                # For non-imaging data or fallback, just clamp the coordinates
                xCoord = clamp(xCoord, 1, imgWidth)
                yCoord = yCoord < 0 ? yCoord : -clamp(yCoord, 1, imgHeight)
            end
            
            selectedTab = "tab2"
            fTime = time()
            eTime = round(fTime - sTime, digits=3)
            msg = "Plot loaded in $(eTime) seconds"
            log_memory_usage("XY Plot Generated", msi_data)
        catch e
            msg = "Could not retrieve spectrum: $e"
            warning_msg = true
            @error "Spectrum plotting failed" exception=(e, catch_backtrace())
        finally
            progressSpectraPlot = false
            btnPlotDisable = false
            btnSpectraDisable = false
            btnStartDisable = false
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end

    # --- Main View Handlers ---
    @onbutton imgMinus begin
        if isempty(selected_folder_main) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_main)
        # Check if folder exists to prevent errors
        if !isdir(folder_path) return end
        
        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=decrement_image(current_msi, msi_bmp)
        new_col_msi=decrement_image(current_col_msi, col_msi_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_msi = new_msi
            current_col_msi = new_col_msi
            imgInt = "/$(selected_folder_main)/$(current_msi)?t=$(timestamp)"
            colorbar = "/$(selected_folder_main)/$(current_col_msi)?t=$(timestamp)"
            text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
            msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImg, plotlayoutImg, _, _ = loadImgPlot(imgInt)
            btnOpticalDisable = false
        end
    end
    @onbutton imgPlus begin
        if isempty(selected_folder_main) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_main)
        if !isdir(folder_path) return end

        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)
        
        new_msi=increment_image(current_msi, msi_bmp)
        new_col_msi=increment_image(current_col_msi, col_msi_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_msi = new_msi
            current_col_msi = new_col_msi
            imgInt = "/$(selected_folder_main)/$(current_msi)?t=$(timestamp)"
            colorbar = "/$(selected_folder_main)/$(current_col_msi)?t=$(timestamp)"
            text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
            msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImg, plotlayoutImg, _, _ = loadImgPlot(imgInt)
            btnOpticalDisable = false
        end
    end

    @onbutton imgMinusT begin
        if isempty(selected_folder_main) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_main)
        if !isdir(folder_path) return end

        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=decrement_image(current_triq, triq_bmp)
        new_col_msi=decrement_image(current_col_triq, col_triq_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_triq = new_msi
            current_col_triq = new_col_msi
            imgIntT = "/$(selected_folder_main)/$(current_triq)?t=$(timestamp)"
            colorbarT = "/$(selected_folder_main)/$(current_col_triq)?t=$(timestamp)"
            text_nmass = replace(current_triq, r"TrIQ_|.bmp" => "")
            msgtriq = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgT, plotlayoutImgT, _, _ = loadImgPlot(imgIntT)
            btnOpticalDisable = false
        end
    end
    @onbutton imgPlusT begin
        if isempty(selected_folder_main) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_main)
        if !isdir(folder_path) return end

        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=increment_image(current_triq, triq_bmp)
        new_col_msi=increment_image(current_col_triq, col_triq_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_triq = new_msi
            current_col_triq = new_col_msi
            imgIntT = "/$(selected_folder_main)/$(current_triq)?t=$(timestamp)"
            colorbarT = "/$(selected_folder_main)/$(current_col_triq)?t=$(timestamp)"
            text_nmass = replace(current_triq, r"TrIQ_|.bmp" => "")
            msgtriq = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgT, plotlayoutImgT, _, _ = loadImgPlot(imgIntT)
            btnOpticalDisable = false
        end
    end

    # --- Compare View Handlers ---
    @onbutton imgMinusCompLeft begin
        if isempty(selected_folder_compare_left) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_left)
        if !isdir(folder_path) return end

        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=decrement_image(current_msiCompLeft, msi_bmp)
        new_col_msi=decrement_image(current_col_msiCompLeft, col_msi_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_msiCompLeft = new_msi
            current_col_msiCompLeft = new_col_msi
            imgIntCompLeft = "/$(selected_folder_compare_left)/$(current_msiCompLeft)?t=$(timestamp)"
            colorbarCompLeft = "/$(selected_folder_compare_left)/$(current_col_msiCompLeft)?t=$(timestamp)"
            text_nmass = replace(current_msiCompLeft, r"MSI_|.bmp" => "")
            msgimgCompLeft = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgCompLeft, plotlayoutImgCompLeft, _, _ = loadImgPlot(imgIntCompLeft)
        end
    end

    @onbutton imgPlusCompLeft begin
        if isempty(selected_folder_compare_left) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_left)
        if !isdir(folder_path) return end

        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)
        
        new_msi=increment_image(current_msiCompLeft, msi_bmp)
        new_col_msi=increment_image(current_col_msiCompLeft, col_msi_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_msiCompLeft = new_msi
            current_col_msiCompLeft = new_col_msi
            imgIntCompLeft = "/$(selected_folder_compare_left)/$(current_msiCompLeft)?t=$(timestamp)"
            colorbarCompLeft = "/$(selected_folder_compare_left)/$(current_col_msiCompLeft)?t=$(timestamp)"
            text_nmass = replace(current_msiCompLeft, r"MSI_|.bmp" => "")
            msgimgCompLeft = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgCompLeft, plotlayoutImgCompLeft, _, _ = loadImgPlot(imgIntCompLeft)
        end
    end

    @onbutton imgMinusTCompLeft begin
        if isempty(selected_folder_compare_left) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_left)
        if !isdir(folder_path) return end

        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=decrement_image(current_triqCompLeft, triq_bmp)
        new_col_msi=decrement_image(current_col_triqCompLeft, col_triq_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_triqCompLeft = new_msi
            current_col_triqCompLeft = new_col_msi
            imgIntTCompLeft = "/$(selected_folder_compare_left)/$(current_triqCompLeft)?t=$(timestamp)"
            colorbarTCompLeft = "/$(selected_folder_compare_left)/$(current_col_triqCompLeft)?t=$(timestamp)"
            text_nmass = replace(current_triqCompLeft, r"TrIQ_|.bmp" => "")
            msgtriqCompLeft = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgTCompLeft, plotlayoutImgTCompLeft, _, _ = loadImgPlot(imgIntTCompLeft)
        end
    end

    @onbutton imgPlusTCompLeft begin
        if isempty(selected_folder_compare_left) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_left)
        if !isdir(folder_path) return end

        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=increment_image(current_triqCompLeft, triq_bmp)
        new_col_msi=increment_image(current_col_triqCompLeft, col_triq_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_triqCompLeft = new_msi
            current_col_triqCompLeft = new_col_msi
            imgIntTCompLeft = "/$(selected_folder_compare_left)/$(current_triqCompLeft)?t=$(timestamp)"
            colorbarTCompLeft = "/$(selected_folder_compare_left)/$(current_col_triqCompLeft)?t=$(timestamp)"
            text_nmass = replace(current_triqCompLeft, r"TrIQ_|.bmp" => "")
            msgtriqCompLeft = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgTCompLeft, plotlayoutImgTCompLeft, _, _ = loadImgPlot(imgIntTCompLeft)
        end
    end

    @onbutton imgMinusCompRight begin
        if isempty(selected_folder_compare_right) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_right)
        if !isdir(folder_path) return end

        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=decrement_image(current_msiCompRight, msi_bmp)
        new_col_msi=decrement_image(current_col_msiCompRight, col_msi_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_msiCompRight = new_msi
            current_col_msiCompRight = new_col_msi
            imgIntCompRight = "/$(selected_folder_compare_right)/$(current_msiCompRight)?t=$(timestamp)"
            colorbarCompRight = "/$(selected_folder_compare_right)/$(current_col_msiCompRight)?t=$(timestamp)"
            text_nmass = replace(current_msiCompRight, r"MSI_|.bmp" => "")
            msgimgCompRight = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgCompRight, plotlayoutImgCompRight, _, _ = loadImgPlot(imgIntCompRight)
        end
    end

    @onbutton imgPlusCompRight begin
        if isempty(selected_folder_compare_right) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_right)
        if !isdir(folder_path) return end

        msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)
        
        new_msi=increment_image(current_msiCompRight, msi_bmp)
        new_col_msi=increment_image(current_col_msiCompRight, col_msi_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_msiCompRight = new_msi
            current_col_msiCompRight = new_col_msi
            imgIntCompRight = "/$(selected_folder_compare_right)/$(current_msiCompRight)?t=$(timestamp)"
            colorbarCompRight = "/$(selected_folder_compare_right)/$(current_col_msiCompRight)?t=$(timestamp)"
            text_nmass = replace(current_msiCompRight, r"MSI_|.bmp" => "")
            msgimgCompRight = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgCompRight, plotlayoutImgCompRight, _, _ = loadImgPlot(imgIntCompRight)
        end
    end

    @onbutton imgMinusTCompRight begin
        if isempty(selected_folder_compare_right) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_right)
        if !isdir(folder_path) return end

        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=decrement_image(current_triqCompRight, triq_bmp)
        new_col_msi=decrement_image(current_col_triqCompRight, col_triq_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_triqCompRight = new_msi
            current_col_triqCompRight = new_col_msi
            imgIntTCompRight = "/$(selected_folder_compare_right)/$(current_triqCompRight)?t=$(timestamp)"
            colorbarTCompRight = "/$(selected_folder_compare_right)/$(current_col_triqCompRight)?t=$(timestamp)"
            text_nmass = replace(current_triqCompRight, r"TrIQ_|.bmp" => "")
            msgtriqCompRight = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgTCompRight, plotlayoutImgTCompRight, _, _ = loadImgPlot(imgIntTCompRight)
        end
    end

    @onbutton imgPlusTCompRight begin
        if isempty(selected_folder_compare_right) return end
        timestamp=string(time_ns())
        folder_path = joinpath("public", selected_folder_compare_right)
        if !isdir(folder_path) return end

        triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)),lt=natural)
        col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)),lt=natural)

        new_msi=increment_image(current_triqCompRight, triq_bmp)
        new_col_msi=increment_image(current_col_triqCompRight, col_triq_png)
        if new_msi !== nothing && new_col_msi !== nothing
            current_triqCompRight = new_msi
            current_col_triqCompRight = new_col_msi
            imgIntTCompRight = "/$(selected_folder_compare_right)/$(current_triqCompRight)?t=$(timestamp)"
            colorbarTCompRight = "/$(selected_folder_compare_right)/$(current_col_triqCompRight)?t=$(timestamp)"
            text_nmass = replace(current_triqCompRight, r"TrIQ_|.bmp" => "")
            msgtriqCompRight = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
            plotdataImgTCompRight, plotlayoutImgTCompRight, _, _ = loadImgPlot(imgIntTCompRight)
        end
    end

    # This handler will now correctly load the first image from the newly selected folder.
    @onchange selected_folder_main begin
        if msi_data !== nothing
            close(msi_data)
        end
        msi_data = nothing
        log_memory_usage("Folder Changed (msi_data cleared)", msi_data)
        GC.gc()
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0)
        end

        if !isempty(selected_folder_main)
            folder_path = joinpath("public", selected_folder_main)
            if !isdir(folder_path) 
                imgInt = ""
                colorbar = ""
                imgIntT = ""
                colorbarT = ""
                msgimg = "Folder not found."
                msgtriq = "Folder not found."
                plotdataImg = [traceImg]
                plotlayoutImg = layoutImg
                plotdataImgT = [traceImg]
                plotlayoutImgT = layoutImg
                imgWidth, imgHeight = 0, 0
                return
            end

            # Handle normal images
            msi_bmp = sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir(folder_path)), lt=natural)
            col_msi_png = sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir(folder_path)), lt=natural)
            
            if !isempty(msi_bmp)
                current_msi = first(msi_bmp)
                imgInt = "/$(selected_folder_main)/$(current_msi)"
                plotdataImg, plotlayoutImg, w, h = loadImgPlot(imgInt)
                imgWidth, imgHeight = w, h
                text_nmass = replace(current_msi, r"MSI_|.bmp" => "")
                msgimg = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
                if !isempty(col_msi_png)
                    current_col_msi = first(col_msi_png)
                    colorbar = "/$(selected_folder_main)/$(current_col_msi)"
                else
                    colorbar = ""
                end
            else
                imgInt = ""
                colorbar = ""
                msgimg = "No MSI images found in this dataset."
                plotdataImg = [traceImg]
                plotlayoutImg = layoutImg
            end

            # Handle TrIQ images
            triq_bmp = sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir(folder_path)), lt=natural)
            col_triq_png = sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir(folder_path)), lt=natural)

            if !isempty(triq_bmp)
                current_triq = first(triq_bmp)
                imgIntT = "/$(selected_folder_main)/$(current_triq)"
                plotdataImgT, plotlayoutImgT, w, h = loadImgPlot(imgIntT)
                # If no MSI image was loaded, dimensions from TrIQ image are used.
                if isempty(msi_bmp)
                    imgWidth, imgHeight = w, h
                end
                text_nmass = replace(current_triq, r"TrIQ_|.bmp" => "")
                msgtriq = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"
                if !isempty(col_triq_png)
                    current_col_triq = first(col_triq_png)
                    colorbarT = "/$(selected_folder_main)/$(current_col_triq)"
                else
                    colorbarT = ""
                end
            else
                imgIntT = ""
                colorbarT = ""
                msgtriq = "No TrIQ images found in this dataset."
                plotdataImgT = [traceImg]
                plotlayoutImgT = layoutImg
            end

            if isempty(msi_bmp) && isempty(triq_bmp)
                imgWidth, imgHeight = 0, 0
            end
        end
    end

    @onchange selected_folder_compare_left begin
        if !isempty(selected_folder_compare_left)
            timestamp = string(time_ns())
            folder_path = joinpath("public", selected_folder_compare_left)
            
            if !isdir(folder_path)
                imgIntCompLeft, colorbarCompLeft, imgIntTCompLeft, colorbarTCompLeft = "", "", "", ""
                msgimgCompLeft, msgtriqCompLeft = "Folder not found.", "Folder not found."
                return
            end

            # Handle normal images
            msi_bmp = sort(filter(f -> startswith(f, "MSI_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
            col_msi_png = sort(filter(f -> startswith(f, "colorbar_MSI_") && endswith(f, ".png"), readdir(folder_path)), lt=natural)
            
            if !isempty(msi_bmp)
                current_msiCompLeft = first(msi_bmp)
                imgIntCompLeft = "/$(selected_folder_compare_left)/$(current_msiCompLeft)?t=$(timestamp)"
                plotdataImgCompLeft, plotlayoutImgCompLeft, _, _ = loadImgPlot(imgIntCompLeft)
                text_nmass = replace(current_msiCompLeft, r"MSI_|.bmp" => "")
                msgimgCompLeft = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
                
                if !isempty(col_msi_png)
                    current_col_msiCompLeft = first(col_msi_png)
                    colorbarCompLeft = "/$(selected_folder_compare_left)/$(current_col_msiCompLeft)?t=$(timestamp)"
                else
                    colorbarCompLeft = ""
                end
            else
                imgIntCompLeft, colorbarCompLeft, msgimgCompLeft = "", "", "No MSI images."
            end

            # Handle TrIQ images
            triq_bmp = sort(filter(f -> startswith(f, "TrIQ_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
            col_triq_png = sort(filter(f -> startswith(f, "colorbar_TrIQ_") && endswith(f, ".png"), readdir(folder_path)), lt=natural)

            if !isempty(triq_bmp)
                current_triqCompLeft = first(triq_bmp)
                imgIntTCompLeft = "/$(selected_folder_compare_left)/$(current_triqCompLeft)?t=$(timestamp)"
                plotdataImgTCompLeft, plotlayoutImgTCompLeft, _, _ = loadImgPlot(imgIntTCompLeft)
                text_nmass = replace(current_triqCompLeft, r"TrIQ_|.bmp" => "")
                msgtriqCompLeft = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"

                if !isempty(col_triq_png)
                    current_col_triqCompLeft = first(col_triq_png)
                    colorbarTCompLeft = "/$(selected_folder_compare_left)/$(current_col_triqCompLeft)?t=$(timestamp)"
                else
                    colorbarTCompLeft = ""
                end
            else
                imgIntTCompLeft, colorbarTCompLeft, msgtriqCompLeft = "", "", "No TrIQ images."
            end
        end
    end

    @onchange selected_folder_compare_right begin
        if !isempty(selected_folder_compare_right)
            timestamp = string(time_ns())
            folder_path = joinpath("public", selected_folder_compare_right)
            if !isdir(folder_path) 
                imgIntCompRight, colorbarCompRight, imgIntTCompRight, colorbarTCompRight = "", "", "", ""
                msgimgCompRight, msgtriqCompRight = "Folder not found.", "Folder not found."
                return
            end

            # Handle normal images
            msi_bmp = sort(filter(f -> startswith(f, "MSI_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
            col_msi_png = sort(filter(f -> startswith(f, "colorbar_MSI_") && endswith(f, ".png"), readdir(folder_path)), lt=natural)
            
            if !isempty(msi_bmp)
                current_msiCompRight = first(msi_bmp)
                imgIntCompRight = "/$(selected_folder_compare_right)/$(current_msiCompRight)?t=$(timestamp)"
                plotdataImgCompRight, plotlayoutImgCompRight, _, _ = loadImgPlot(imgIntCompRight)
                text_nmass = replace(current_msiCompRight, r"MSI_|.bmp" => "")
                msgimgCompRight = "<i>m/z</i>: $(replace(text_nmass, "_" => "."))"
                
                if !isempty(col_msi_png)
                    current_col_msiCompRight = first(col_msi_png)
                    colorbarCompRight = "/$(selected_folder_compare_right)/$(current_col_msiCompRight)?t=$(timestamp)"
                else
                    colorbarCompRight = ""
                end
            else
                imgIntCompRight, colorbarCompRight, msgimgCompRight = "", "", "No MSI images."
            end

            # Handle TrIQ images
            triq_bmp = sort(filter(f -> startswith(f, "TrIQ_") && endswith(f, ".bmp"), readdir(folder_path)), lt=natural)
            col_triq_png = sort(filter(f -> startswith(f, "colorbar_TrIQ_") && endswith(f, ".png"), readdir(folder_path)), lt=natural)

            if !isempty(triq_bmp)
                current_triqCompRight = first(triq_bmp)
                imgIntTCompRight = "/$(selected_folder_compare_right)/$(current_triqCompRight)?t=$(timestamp)"
                plotdataImgTCompRight, plotlayoutImgTCompRight, _, _ = loadImgPlot(imgIntTCompRight)
                text_nmass = replace(current_triqCompRight, r"TrIQ_|.bmp" => "")
                msgtriqCompRight = "TrIQ <i>m/z</i>: $(replace(text_nmass, "_" => "."))"

                if !isempty(col_triq_png)
                    current_col_triqCompRight = first(col_triq_png)
                    colorbarTCompRight = "/$(selected_folder_compare_right)/$(current_col_triqCompRight)?t=$(timestamp)"
                else
                    colorbarTCompRight = ""
                end
            else
                imgIntTCompRight, colorbarTCompRight, msgtriqCompRight = "", "", "No TrIQ images."
            end
        end
    end

    
    # 3d plot
    @onbutton image3dPlot begin
        msg = "Image 3D plot selected"
        cleaned_imgInt = replace(imgInt, r"\?.*" => "")
        cleaned_imgInt = lstrip(cleaned_imgInt, '/')
        var = joinpath("./public", cleaned_imgInt)

        if !isfile(var)
            msg = "Image could not be 3d plotted"
            warning_msg = true
            return
        end

        progressPlot = true
        btnPlotDisable = true
        btnStartDisable = true
        btnSpectraDisable = true

        try
            # --- Get Mask Path ---
            local mask_path_for_plot::Union{String, Nothing} = nothing
            if maskEnabled && !isempty(selected_folder_main)
                registry = load_registry(registry_path)
                entry = get(registry, selected_folder_main, nothing)
                if entry !== nothing && get(entry, "has_mask", false)
                    mask_path_candidate = get(entry, "mask_path", "")
                    if isfile(mask_path_candidate)
                        mask_path_for_plot = mask_path_candidate
                    else
                        @warn "Mask enabled but file not found: $(mask_path_candidate). Plotting without mask."
                    end
                end
            end
            # ---

            sTime = time()
            if mask_path_for_plot !== nothing
                plotdata3d, plotlayout3d = loadSurfacePlot(imgInt, mask_path_for_plot)
            else
                plotdata3d, plotlayout3d = loadSurfacePlot(imgInt)
            end

            selectedTab = "tab4"
            fTime = time()
            eTime = round(fTime - sTime, digits=3)
            msg = "Plot loaded in $(eTime) seconds"
            log_memory_usage("Mean Plot Generated", msi_data)
        catch e
            msg = "Failed to load and process image: $e"
            warning_msg = true
            @error "3D plot generation failed" exception=(e, catch_backtrace())
        finally
            progressPlot=false
            btnPlotDisable=false
            btnStartDisable=false
            btnSpectraDisable=false
            SpectraEnabled=true
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end
    
    @onbutton triq3dPlot begin
        msg = "TrIQ 3D plot selected"
        cleaned_imgIntT = replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT = lstrip(cleaned_imgIntT, '/')
        var = joinpath("./public", cleaned_imgIntT)

        if !isfile(var)
            msg = "Image could not be 3d plotted"
            warning_msg = true
            return
        end

        progressPlot = true
        btnPlotDisable = true
        btnStartDisable = true
        btnSpectraDisable = true

        try
            # --- Get Mask Path ---
            local mask_path_for_plot::Union{String, Nothing} = nothing
            if maskEnabled[] && !isempty(selected_folder_main)
                registry = load_registry(registry_path)
                entry = get(registry, selected_folder_main, nothing)
                if entry !== nothing && get(entry, "has_mask", false)
                    mask_path_candidate = get(entry, "mask_path", "")
                    if isfile(mask_path_candidate)
                        mask_path_for_plot = mask_path_candidate
                    else
                        @warn "Mask enabled but file not found: $(mask_path_candidate). Plotting without mask."
                    end
                end
            end
            # ---

            sTime = time()
            if mask_path_for_plot !== nothing
                plotdata3d, plotlayout3d = loadSurfacePlot(imgIntT, mask_path_for_plot)
            else
                plotdata3d, plotlayout3d = loadSurfacePlot(imgIntT)
            end

            selectedTab = "tab4"
            fTime = time()
            eTime = round(fTime - sTime, digits=3)
            msg = "Plot loaded in $(eTime) seconds"
            log_memory_usage("Mean Plot Generated", msi_data)
        catch e
            msg = "Failed to load and process image: $e"
            warning_msg = true
            @error "3D TrIQ plot generation failed" exception=(e, catch_backtrace())
        finally
            progressPlot=false
            btnPlotDisable=false
            btnStartDisable=false
            btnSpectraDisable=false
            SpectraEnabled=true
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end

    # Contour 2d plot
    @onbutton imageCPlot begin
        msg="Image 2D plot selected"
        cleaned_imgInt=replace(imgInt, r"\?.*" => "")
        cleaned_imgInt=lstrip(cleaned_imgInt, '/')
        var=joinpath("./public", cleaned_imgInt)
        
        if !isfile(var)
            msg="Image could not be 2D plotted"
            warning_msg=true
            return
        end

        progressPlot=true
        btnPlotDisable=true
        btnStartDisable=true
        btnSpectraDisable=true
        
        try
            sTime=time()
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
        finally
            progressPlot=false
            btnPlotDisable=false
            btnStartDisable=false
            btnSpectraDisable=false
            SpectraEnabled=true
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end
    # Contour 2d plot for TrIQ
    @onbutton triqCPlot begin
        msg="Image 2D plot selected"
        cleaned_imgIntT=replace(imgIntT, r"\?.*" => "")
        cleaned_imgIntT=lstrip(cleaned_imgIntT, '/')
        var=joinpath("./public", cleaned_imgIntT)
        
        if !isfile(var)
            msg="Image could not be 2D plotted"
            warning_msg=true
            return
        end

        progressPlot=true
        btnPlotDisable=true
        btnStartDisable=true
        btnSpectraDisable=true
        
        try
            sTime=time()
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
        finally
            progressPlot=false
            btnPlotDisable=false
            btnStartDisable=false
            btnSpectraDisable=false
            SpectraEnabled=true
            GC.gc()
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0)
            end
        end
    end

    @onbutton compareBtn begin
        CompareDialog=true
    end

    # To include a visualization in the spectrum plot indicating where is the selected mass
    @onchange Nmass begin
        if !isempty(xSpectraMz)
            # Main spectrum trace
            traceSpectra = PlotlyBase.stem(
                x=xSpectraMz, 
                y=ySpectraMz, 
                marker=attr(size=1, color="blue", opacity=0.5), 
                name="Spectrum", 
                hoverinfo="x", 
                hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>", 
                showlegend=false
            )
            
            # Parse all valid masses from the comma-separated string
            mass_strs = split(Nmass, ',', keepempty=false)
            mass_traces = [traceSpectra]  # Start with the main spectrum
            
            valid_masses = Float64[]
            for (idx, mass_str) in enumerate(mass_strs)
                try
                    mass_val = parse(Float64, strip(mass_str))
                    if mass_val > 0  # Only add valid positive masses
                        push!(valid_masses, mass_val)
                        
                        # Create a vertical line for this mass (Plotly will auto-assign colors)
                        mass_trace = PlotlyBase.scatter(
                            x=[mass_val, mass_val], 
                            y=[0, maximum(ySpectraMz)], 
                            mode="lines", 
                            line=attr(width=1.5, dash="dash"),
                            name="m/z $(round(mass_val, digits=4))",
                            showlegend=false,
                            hoverinfo="x+name",
                            hovertemplate="<b>%{data.name}</b><extra></extra>"
                        )
                        push!(mass_traces, mass_trace)
                    end
                catch e
                    # Skip invalid entries, continue with next
                    continue
                end
            end
            
            # Update the plot data
            plotdata = mass_traces
        end
    end

    # Event detection for clicking on the images
    @onchange data_click begin
        if selectedTab == "tab1" || selectedTab == "tab0"
            # This is for the image heatmaps
            cursor_data = get(data_click, "cursor", nothing)
            if cursor_data === nothing
                return
            end

            x_val = get(cursor_data, "x", nothing)
            y_val = get(cursor_data, "y", nothing)

            if x_val === nothing || y_val === nothing
                return # Do nothing if coordinates are not provided by the event
            end

            x = Int32(round(x_val))
            y = Int32(round(y_val)) # y is negative in the UI

            # Update the reactive coordinates, which will trigger the crosshair update
            xCoord = clamp(x, 1, imgWidth)
            yCoord = clamp(y, -imgHeight, -1)
        end
    end

    @onchange xCoord, yCoord begin 
        if selectedTab == "tab1"
            main_trace = plotdataImgT[1]  # The heatmap/image trace
            trace1, trace2 = crossLinesPlot(xCoord, yCoord, imgWidth, -imgHeight)
            plotdataImgT = [main_trace, trace1, trace2]  # Fresh array every time
        elseif selectedTab == "tab0"
            main_trace = plotdataImg[1]  # The heatmap/image trace
            trace1, trace2 = crossLinesPlot(xCoord, yCoord, imgWidth, -imgHeight)
            plotdataImg = [main_trace, trace1, trace2]
        end
    end

    @onbutton btnOptical begin
        imgRoute=pick_file(; filterlist="png,bmp,jpg,jpeg")
        if imgRoute==""
            msg="No optical image selected"
        else
            selectedTab="tab0"
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
            img=load(imgRoute)
            save("./public/css/imgOver.png",img)
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt,"/css/imgOver.png",imgTrans)
        end
        
    end

    @onbutton btnOpticalT begin
        imgRoute=pick_file(; filterlist="png,bmp,jpg,jpeg")
        if imgRoute==""
            msg="No optical image selected"
        else
            selectedTab="tab1"
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
            img=load(imgRoute)
            save("./public/css/imgOver.png",img)
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT,"/css/imgOver.png",imgTrans)
            opticalOverTriq=true
        end
    end

    @onchange imgTrans begin
        if !opticalOverTriq && imgRoute!=""
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt,"/css/imgOver.png",imgTrans)
        elseif opticalOverTriq && imgRoute!=""
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT,"/css/imgOver.png",imgTrans)
        end
    end

    @onchange opticalOverTriq begin
        if !opticalOverTriq && imgRoute!=""
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt,"/css/imgOver.png",imgTrans)
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT)
            selectedTab="tab0"
        elseif opticalOverTriq && imgRoute!=""
            plotdataImg, plotlayoutImg, imgWidth, imgHeight=loadImgPlot(imgInt)
            plotdataImgT, plotlayoutImgT, imgWidth, imgHeight=loadImgPlot(imgIntT,"/css/imgOver.png",imgTrans)
            selectedTab="tab1"
        end
    end

    @mounted watchplots()

    @onchange isready @time begin
        if isready && !registry_init_done
            warmup_init()
            try
                println("Synchronizing registry with filesystem on backend init...")
                reg_path = abspath(joinpath(@__DIR__, "public", "registry.json"))
                registry = isfile(reg_path) ? load_registry(reg_path) : Dict{String, Any}()

                public_dirs = isdir("public") ? readdir("public") : []
                ignored_dirs = ["css", "masks"]
                
                dataset_dirs = filter(d -> isdir(joinpath("public", d)) && !(d in ignored_dirs), public_dirs)
                
                registry_keys = Set(keys(registry))
                folder_set = Set(dataset_dirs)

                new_folders = setdiff(folder_set, registry_keys)
                for folder in new_folders
                    println("Found new folder: $folder")
                    registry[folder] = Dict(
                        "source_path" => "unknown (manually added)",
                        "processed_date" => "unknown",
                        "metadata" => Dict(),
                        "is_imzML" => true # Assume folder contains images if found this way
                    )
                end

                removed_folders = setdiff(registry_keys, folder_set)
                for folder in removed_folders
                    delete!(registry, folder)
                end

                if !isempty(new_folders) || !isempty(removed_folders)
                    println("Registry changed, saving...")
                    open(reg_path, "w") do f
                        JSON.print(f, registry, 4)
                    end
                end
                
                all_folders = sort(collect(keys(registry)), lt=natural)
                img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)

                available_folders = deepcopy(all_folders)
                image_available_folders = deepcopy(img_folders)
                
                println("UI lists updated. All: $(length(available_folders)), Images: $(length(image_available_folders))")

            catch e
                @warn "Registry synchronization failed: $e"
                available_folders = []
                image_available_folders = []
                selected_files = String[]
            finally
                registry_init_done = true
            end
        end
        log_memory_usage("App Ready", msi_data)
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
