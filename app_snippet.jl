# app.jl

module App
# ==Packages ==
using GenieFramework
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
using .MSI_src: MSIData, OpenMSIData, process_spectrum, IterateSpectra, ImzMLSource, _iterate_spectra_fast, MzMLSource, find_mass, ViridisPalette, get_mz_slice, get_multiple_mz_slices, quantize_intensity, save_bitmap, median_filter, save_bitmap, downsample_spectrum, TrIQ, precompute_analytics, ImportMzmlFile, generate_colorbar_image, load_and_prepare_mask, set_global_mz_range!, main_precalculation, MutableSpectrum

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
    # == Loading Screen Variables ==
    @in is_initializing = true
    @in initialization_message = "Initializing..."

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

    ## Preprocessing Parameters
    @in progressPrep=false
    @in stabilization_method="sqrt"
    @in smoothing_method="sg"
    @in smoothing_window = ""
    @in smoothing_order = ""
    @in baseline_method="snip"
    @in baseline_iterations = ""
    @in baseline_window = ""
    @in normalization_method="tic"
    @in alignment_method="lowess"
    @in alignment_span = ""
    @in alignment_tolerance = ""
    @in alignment_tolerance_unit="mz"
    @in alignment_max_shift_ppm = ""
    @in alignment_min_matched_peaks = ""
    @in peak_picking_method="profile"
    @in peak_picking_snr_threshold = ""
    @in peak_picking_half_window = ""
    @in peak_picking_min_peak_prominence = ""
    @in peak_picking_merge_peaks_tolerance = ""
    @in peak_picking_min_peak_width_ppm = ""
    @in peak_picking_max_peak_width_ppm = ""
    @in peak_picking_min_peak_shape_r2 = ""
    @in binning_method="adaptive"
    @in binning_tolerance = ""
    @in binning_tolerance_unit="ppm"
    @in binning_frequency_threshold = ""
    @in binning_min_peak_per_bin = ""
    @in binning_max_bin_width_ppm = ""
    @in binning_intensity_weighted_centers=true
    @in binning_num_uniform_bins = ""
    @in calibration_fit_order = ""
    @in calibration_ppm_tolerance = ""
    @in peak_selection_min_snr = ""
    @in peak_selection_min_fwhm_ppm = ""
    @in peak_selection_max_fwhm_ppm = ""
    @in peak_selection_min_shape_r2 = ""
    @in peak_selection_frequency_threshold = ""
    @in peak_selection_correlation_threshold = ""

    @in reference_peaks_list = [
        Dict("mz" => 137.0244, "label" => "DHB_fragment"),
        Dict("mz" => 155.0349, "label" => "DHB_M+H"),
    ]

    # --- Methods for Reference Peaks List ---
    function addReferencePeak()
        push!(reference_peaks_list, Dict("mz" => 0.0, "label" => ""))
        reference_peaks_list = deepcopy(reference_peaks_list) # Force reactivity
    end

    function removeReferencePeak(index::Int)
        deleteat!(reference_peaks_list, index)
        reference_peaks_list = deepcopy(reference_peaks_list) # Force reactivity
    end

    # Individual step enable/disable flags
    @in enable_stabilization = true
    @in enable_smoothing = true
    @in enable_baseline = true
    @in enable_normalization = true
    @in enable_standards = true
    @in enable_alignment = true
    @in enable_peak_picking = true
    @in enable_binning = true
    @in enable_calibration = true
    @in enable_peak_selection = true

    # Trigger for running the full pipeline
    @in run_full_pipeline = false
    @out current_pipeline_step = "" # To indicate which step is currently running in the full pipeline

    @in export_params_btn = false
    @in import_params_btn = false
    @in imported_params_file = nothing

    @in suggested_smoothing_window = ""
    @in suggested_smoothing_order = ""
    @in suggested_baseline_iterations = ""
    @in suggested_baseline_window = ""
    @in suggested_alignment_span = ""
    @in suggested_alignment_tolerance = ""
    @in suggested_alignment_max_shift_ppm = ""
    @in suggested_alignment_min_matched_peaks = ""
    @in suggested_peak_picking_snr_threshold = ""
    @in suggested_peak_picking_half_window = ""
    @in suggested_peak_picking_min_peak_prominence = ""
    @in suggested_peak_picking_merge_peaks_tolerance = ""
    @in suggested_peak_picking_min_peak_width_ppm = ""
    @in suggested_peak_picking_max_peak_width_ppm = ""
    @in suggested_peak_picking_min_peak_shape_r2 = ""
    @in suggested_binning_tolerance = ""
    @in suggested_binning_frequency_threshold = ""
    @in suggested_binning_min_peak_per_bin = ""
    @in suggested_binning_max_bin_width_ppm = ""
    @in suggested_binning_num_uniform_bins = ""
    @in suggested_calibration_fit_order = ""
    @in suggested_calibration_ppm_tolerance = ""
    @in suggested_peak_selection_min_snr = ""
    @in suggested_peak_selection_min_fwhm_ppm = ""
    @in suggested_peak_selection_max_fwhm_ppm = ""
    @in suggested_peak_selection_min_shape_r2 = ""
    @in suggested_peak_selection_frequency_threshold = ""
    @in suggested_peak_selection_correlation_threshold = ""

    # == Batch Summary Dialog ==
    @in showBatchSummary = false
    @out batch_summary = ""

    # == Batch Processing & Registry Variables ==
    @private registry_init_done = false
    @in refetch_folders = false
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

    # UI plot data for Preprocessing
    @out plotdata_before = [traceSpectra]
    @out plotlayout_before = layoutSpectra
    @out plotdata_after = [traceSpectra]
    @out plotlayout_after = layoutSpectra

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
            dataset_name = replace(basename(picked_route), r"(\.(imzML|imzml|mzML|mzml))$ "i => "")
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
            
            # Auto-suggest parameters
            try
                println("Calling main_precalculation to get recommended parameters...")
                recommended_params = main_precalculation(loaded_data)
                
                for (step_name, params) in recommended_params
                    for (param_key, value) in params
                        # Convert value to appropriate type before assignment
                        processed_value = if value === nothing
                            nothing
                        elseif value isa Tuple
                            @warn "Skipping invalid parameter suggestion (tuple): $value for $param_key"
                            "" # Set to empty string for safety
                        elseif value isa Number
                            value
                        else
                            string(value)
                        end
                        
                        if processed_value !== nothing
                            if step_name == :Smoothing
                                if param_key == :window
                                    suggested_smoothing_window = string(processed_value)
                                    smoothing_window = string(processed_value)
                                    println("    suggested_smoothing_window set to $(suggested_smoothing_window)")
                                elseif param_key == :order
                                    suggested_smoothing_order = string(processed_value)
                                    smoothing_order = string(processed_value)
                                    println("    suggested_smoothing_order set to $(suggested_smoothing_order)")
                                end
                            elseif step_name == :BaselineCorrection
                                if param_key == :iterations
                                    suggested_baseline_iterations = string(processed_value)
                                    baseline_iterations = string(processed_value)
                                    println("    suggested_baseline_iterations set to $(suggested_baseline_iterations)")
                                elseif param_key == :window
                                    suggested_baseline_window = string(processed_value)
                                    baseline_window = string(processed_value)
                                    println("    suggested_baseline_window set to $(suggested_baseline_window)")
                                end
                            elseif step_name == :PeakAlignment
                                if param_key == :span
                                    suggested_alignment_span = string(processed_value)
                                    alignment_span = string(processed_value)
                                    println("    suggested_alignment_span set to $(suggested_alignment_span)")
                                elseif param_key == :tolerance
                                    suggested_alignment_tolerance = string(processed_value)
                                    alignment_tolerance = string(processed_value)
                                    println("    suggested_alignment_tolerance set to $(suggested_alignment_tolerance)")
                                elseif param_key == :max_shift_ppm
                                    suggested_alignment_max_shift_ppm = string(processed_value)
                                    alignment_max_shift_ppm = string(processed_value)
                                    println("    suggested_alignment_max_shift_ppm set to $(suggested_alignment_max_shift_ppm)")
                                elseif param_key == :min_matched_peaks
                                    suggested_alignment_min_matched_peaks = string(processed_value)
                                    alignment_min_matched_peaks = string(processed_value)
                                    println("    suggested_alignment_min_matched_peaks set to $(suggested_alignment_min_matched_peaks)")
                                end
                            elseif step_name == :Calibration
                                if param_key == :fit_order
                                    suggested_calibration_fit_order = string(processed_value)
                                    calibration_fit_order = string(processed_value)
                                    println("    suggested_calibration_fit_order set to $(suggested_calibration_fit_order)")
                                elseif param_key == :ppm_tolerance
                                    suggested_calibration_ppm_tolerance = string(processed_value)
                                    calibration_ppm_tolerance = string(processed_value)
                                    println("    suggested_calibration_ppm_tolerance set to $(suggested_calibration_ppm_tolerance)")
                                end
                            elseif step_name == :PeakPicking
                                if param_key == :snr_threshold
                                    suggested_peak_picking_snr_threshold = string(processed_value)
                                    peak_picking_snr_threshold = string(processed_value)
                                    println("    suggested_peak_picking_snr_threshold set to $(suggested_peak_picking_snr_threshold)")
                                elseif param_key == :half_window
                                    suggested_peak_picking_half_window = string(processed_value)
                                    peak_picking_half_window = string(processed_value)
                                    println("    suggested_peak_picking_half_window set to $(suggested_peak_picking_half_window)")
                                elseif param_key == :min_peak_prominence
                                    suggested_peak_picking_min_peak_prominence = string(processed_value)
                                    peak_picking_min_peak_prominence = string(processed_value)
                                    println("    suggested_peak_picking_min_peak_prominence set to $(suggested_peak_picking_min_peak_prominence)")
                                elseif param_key == :merge_peaks_tolerance
                                    suggested_peak_picking_merge_peaks_tolerance = string(processed_value)
                                    peak_picking_merge_peaks_tolerance = string(processed_value)
                                    println("    suggested_peak_picking_merge_peaks_tolerance set to $(suggested_peak_picking_merge_peaks_tolerance)")
                                elseif param_key == :min_peak_width_ppm
                                    suggested_peak_picking_min_peak_width_ppm = string(processed_value)
                                    peak_picking_min_peak_width_ppm = string(processed_value)
                                    println("    suggested_peak_picking_min_peak_width_ppm set to $(suggested_peak_picking_min_peak_width_ppm)")
                                elseif param_key == :max_peak_width_ppm
                                    suggested_peak_picking_max_peak_width_ppm = string(processed_value)
                                    peak_picking_max_peak_width_ppm = string(processed_value)
                                    println("    suggested_peak_picking_max_peak_width_ppm set to $(suggested_peak_picking_max_peak_width_ppm)")
                                elseif param_key == :min_peak_shape_r2
                                    suggested_peak_picking_min_peak_shape_r2 = string(processed_value)
                                    peak_picking_min_peak_shape_r2 = string(processed_value)
                                    println("    suggested_peak_picking_min_peak_shape_r2 set to $(suggested_peak_picking_min_peak_shape_r2)")
                                end
                            elseif step_name == :PeakSelection
                                if param_key == :min_snr
                                    suggested_peak_selection_min_snr = string(processed_value)
                                    peak_selection_min_snr = string(processed_value)
                                    println("    suggested_peak_selection_min_snr set to $(suggested_peak_selection_min_snr)")
                                elseif param_key == :min_fwhm_ppm
                                    suggested_peak_selection_min_fwhm_ppm = string(processed_value)
                                    peak_selection_min_fwhm_ppm = string(processed_value)
                                    println("    suggested_peak_selection_min_fwhm_ppm set to $(suggested_peak_selection_min_fwhm_ppm)")
                                elseif param_key == :max_fwhm_ppm
                                    suggested_peak_selection_max_fwhm_ppm = string(processed_value)
                                    peak_selection_max_fwhm_ppm = string(processed_value)
                                    println("    suggested_peak_selection_max_fwhm_ppm set to $(suggested_peak_selection_max_fwhm_ppm)")
                                elseif param_key == :min_shape_r2
                                    suggested_peak_selection_min_shape_r2 = string(processed_value)
                                    peak_selection_min_shape_r2 = string(processed_value)
                                    println("    suggested_peak_selection_min_shape_r2 set to $(suggested_peak_selection_min_shape_r2)")
                                elseif param_key == :frequency_threshold
                                    suggested_peak_selection_frequency_threshold = string(processed_value)
                                    peak_selection_frequency_threshold = string(processed_value)
                                    println("    suggested_peak_selection_frequency_threshold set to $(suggested_peak_selection_frequency_threshold)")
                                elseif param_key == :correlation_threshold
                                    suggested_peak_selection_correlation_threshold = string(processed_value)
                                    peak_selection_correlation_threshold = string(processed_value)
                                    println("    suggested_peak_selection_correlation_threshold set to $(suggested_peak_selection_correlation_threshold)")
                                end
                            elseif step_name == :PeakBinning
                                if param_key == :tolerance
                                    suggested_binning_tolerance = string(processed_value)
                                    binning_tolerance = string(processed_value)
                                    println("    suggested_binning_tolerance set to $(suggested_binning_tolerance)")
                                elseif param_key == :frequency_threshold
                                    suggested_binning_frequency_threshold = string(processed_value)
                                    binning_frequency_threshold = string(processed_value)
                                    println("    suggested_binning_frequency_threshold set to $(suggested_binning_frequency_threshold)")
                                elseif param_key == :min_peak_per_bin
                                    suggested_binning_min_peak_per_bin = string(processed_value)
                                    binning_min_peak_per_bin = string(processed_value)
                                    println("    suggested_binning_min_peak_per_bin set to $(suggested_binning_min_peak_per_bin)")
                                elseif param_key == :max_bin_width_ppm
                                    suggested_binning_max_bin_width_ppm = string(processed_value)
                                    binning_max_bin_width_ppm = string(processed_value)
                                    println("    suggested_binning_max_bin_width_ppm set to $(suggested_binning_max_bin_width_ppm)")
                                elseif param_key == :num_uniform_bins
                                    suggested_binning_num_uniform_bins = string(processed_value)
                                    binning_num_uniform_bins = string(processed_value)
                                    println("    suggested_binning_num_uniform_bins set to $(suggested_binning_num_uniform_bins)")
                                end
                            end
                        end
                    end
                end

                # Also set method types for steps
                if haskey(recommended_params, :Smoothing) && haskey(recommended_params[:Smoothing], :method)
                    smoothing_method = string(recommended_params[:Smoothing][:method])
                end
                if haskey(recommended_params, :BaselineCorrection) && haskey(recommended_params[:BaselineCorrection], :method)
                    baseline_method = string(recommended_params[:BaselineCorrection][:method])
                end
                if haskey(recommended_params, :Normalization) && haskey(recommended_params[:Normalization], :method)
                    normalization_method = string(recommended_params[:Normalization][:method])
                end
                if haskey(recommended_params, :PeakAlignment) && haskey(recommended_params[:PeakAlignment], :method)
                    alignment_method = string(recommended_params[:PeakAlignment][:method])
                end
                if haskey(recommended_params, :PeakPicking) && haskey(recommended_params[:PeakPicking], :method)
                    peak_picking_method = string(recommended_params[:PeakPicking][:method])
                end
                if haskey(recommended_params, :PeakBinning) && haskey(recommended_params[:PeakBinning], :method)
                    binning_method = string(recommended_params[:PeakBinning][:method])
                end
                if haskey(recommended_params, :PeakAlignment) && haskey(recommended_params[:PeakAlignment], :tolerance_unit)
                    alignment_tolerance_unit = string(recommended_params[:PeakAlignment][:tolerance_unit])
                end
                if haskey(recommended_params, :PeakBinning) && haskey(recommended_params[:PeakBinning], :tolerance_unit)
                    binning_tolerance_unit = string(recommended_params[:PeakBinning][:tolerance_unit])
                end

                msg = "File loaded and parameters suggested."
            catch e
                @warn "Could not suggest parameters. Using defaults. Error: $e"
            end

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
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
            end
            progress = false
            progressSpectraPlot = false
        end
    end

    @onbutton export_params_btn begin
        params_to_export = Dict(
            "stabilization_method" => stabilization_method,
            "smoothing_method" => smoothing_method,
            "smoothing_window" => smoothing_window,
            "smoothing_order" => smoothing_order,
            "baseline_method" => baseline_method,
            "baseline_iterations" => baseline_iterations,
            "baseline_window" => baseline_window,
            "normalization_method" => normalization_method,
            "alignment_method" => alignment_method,
            "alignment_span" => alignment_span,
            "alignment_tolerance" => alignment_tolerance,
            "alignment_tolerance_unit" => alignment_tolerance_unit,
            "alignment_max_shift_ppm" => alignment_max_shift_ppm,
            "alignment_min_matched_peaks" => alignment_min_matched_peaks,
            "peak_picking_method" => peak_picking_method,
            "peak_picking_snr_threshold" => peak_picking_snr_threshold,
            "peak_picking_half_window" => peak_picking_half_window,
            "peak_picking_min_peak_prominence" => peak_picking_min_peak_prominence,
            "peak_picking_merge_peaks_tolerance" => peak_picking_merge_peaks_tolerance,
            "peak_picking_min_peak_width_ppm" => peak_picking_min_peak_width_ppm,
            "peak_picking_max_peak_width_ppm" => peak_picking_max_peak_width_ppm,
            "peak_picking_min_peak_shape_r2" => peak_picking_min_peak_shape_r2,
            "binning_method" => binning_method,
            "binning_tolerance" => binning_tolerance,
            "binning_tolerance_unit" => binning_tolerance_unit,
            "binning_frequency_threshold" => binning_frequency_threshold,
            "binning_min_peak_per_bin" => binning_min_peak_per_bin,
            "binning_max_bin_width_ppm" => binning_max_bin_width_ppm,
            "binning_intensity_weighted_centers" => binning_intensity_weighted_centers,
            "binning_num_uniform_bins" => binning_num_uniform_bins,
            "calibration_fit_order" => calibration_fit_order,
            "calibration_ppm_tolerance" => calibration_ppm_tolerance,
            "peak_selection_min_snr" => peak_selection_min_snr,
            "peak_selection_min_fwhm_ppm" => peak_selection_min_fwhm_ppm,
            "peak_selection_max_fwhm_ppm" => peak_selection_max_fwhm_ppm,
            "peak_selection_min_shape_r2" => peak_selection_min_shape_r2,
            "peak_selection_frequency_threshold" => peak_selection_frequency_threshold,
            "peak_selection_correlation_threshold" => peak_selection_correlation_threshold,
            "reference_peaks_list" => reference_peaks_list
        )
        json_string = JSON.json(params_to_export)
        js_script = """
            var element = document.createElement('a');
            element.setAttribute('href', 'data:text/json;charset=utf-8,' + encodeURIComponent(`$json_string`));
            element.setAttribute('download', 'preprocessing_params.json');
            element.style.display = 'none';
            document.body.appendChild(element);
            element.click();
            document.body.removeChild(element);
        """
        run_js(js_script)
        msg = "Parameters exported."
    end

    @onchange import_params_btn begin
        if import_params_btn
            try
                json_string = String(imported_params_file.data)
                params = JSON.parse(json_string)

                for (key, value) in params
                    if key == "reference_peaks_list"
                        reference_peaks_list = value
                    else
                        # Use getfield and setproperty! to update reactive variables by name
                        if hasfield(typeof(@__MODULE__), Symbol(key))
                            getfield(@__MODULE__, Symbol(key))[] = value
                        end
                    end
                end
                msg = "Parameters imported successfully."
            catch e
                msg = "Failed to import parameters: $e"
                warning_msg = true
            finally
                import_params_btn = false # Reset button state
            end
        end
    end

    @onbutton run_full_pipeline begin
        progressPrep = true
        msg = "Running preprocessing pipeline..."
        try
            # 1. Get data from "before" plot
            if isempty(plotdata_before.traces) || isempty(plotdata_before.traces[1].x)
                msg = "Please generate a spectrum plot first (e.g., Mean, Sum, or X,Y)."
                warning_msg = true
                return
            end
            
            mz_data = plotdata_before.traces[1].x
            intensity_data = plotdata_before.traces[1].y
            
            # 2. Create a temporary MutableSpectrum
            temp_spectrum = MutableSpectrum(mz_data, intensity_data, [], 1) # id=1 is arbitrary
            
            # 3. Run the pipeline
            pipeline_spectra = [temp_spectrum]

            if enable_stabilization
                current_pipeline_step = "Stabilizing..."
                params = Dict(:method => Symbol(stabilization_method))
                apply_intensity_transformation(pipeline_spectra, params)
            end
            if enable_smoothing
                current_pipeline_step = "Smoothing..."
                params = Dict(:method => Symbol(smoothing_method), :window => parse(Int, smoothing_window), :order => parse(Int, smoothing_order))
                apply_smoothing(pipeline_spectra, params)
            end
            if enable_baseline
                current_pipeline_step = "Correcting baseline..."
                params = Dict(:method => Symbol(baseline_method), :iterations => parse(Int, baseline_iterations), :window => parse(Int, baseline_window))
                apply_baseline_correction(pipeline_spectra, params)
            end
            if enable_normalization
                current_pipeline_step = "Normalizing..."
                params = Dict(:method => Symbol(normalization_method))
                apply_normalization(pipeline_spectra, params)
            end
            if enable_standards && enable_calibration
                current_pipeline_step = "Calibrating..."
                ref_peaks = Dict(p["mz"] => p["label"] for p in reference_peaks_list)
                params = Dict(:ppm_tolerance => parse(Float64, calibration_ppm_tolerance), :fit_order => parse(Int, calibration_fit_order))
                apply_calibration(pipeline_spectra, params)
            end
            if enable_alignment
                 # Alignment is a no-op for a single spectrum, but we call it for completeness
                current_pipeline_step = "Aligning (skipped for single spectrum)..."
            end
            if enable_peak_picking
                current_pipeline_step = "Picking peaks..."
                params = Dict(
                    :method => Symbol(peak_picking_method), 
                    :snr_threshold => parse(Float64, peak_picking_snr_threshold),
                    :half_window => parse(Int, peak_picking_half_window),
                    :min_peak_prominence => parse(Float64, peak_picking_min_peak_prominence),
                    :merge_peaks_tolerance => parse(Float64, peak_picking_merge_peaks_tolerance)
                )
                apply_peak_picking(pipeline_spectra, params)
            end

            # 4. Update "After" plot
            processed_spectrum = pipeline_spectra[1]
            
            # Main spectrum trace
            after_trace = PlotlyBase.scatter(x=processed_spectrum.mz, y=processed_spectrum.intensity, mode="lines", name="Processed Spectrum")
            
            traces_after = [after_trace]

            # Add peaks if they exist
            if !isempty(processed_spectrum.peaks)
                peak_mzs = [p.mz for p in processed_spectrum.peaks]
                peak_intensities = [p.intensity for p in processed_spectrum.peaks]
                peak_trace = PlotlyBase.scatter(x=peak_mzs, y=peak_intensities, mode="markers", name="Picked Peaks", marker=attr(color="red", size=8))
                push!(traces_after, peak_trace)
            end
            
            plotdata_after = traces_after
            plotlayout_after = PlotlyBase.Layout(title="After Preprocessing")
            msg = "Pipeline finished."

        catch e
            msg = "Error during pipeline execution: $e"
            warning_msg = true
            @error "Pipeline failed" exception=(e, catch_backtrace())
        finally
            progressPrep = false
            current_pipeline_step = ""
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
                progress_message = "Processing file $(file_idx)/$(num_files): $(basename(file_path))"
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
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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
            plotdata_before = plotdata
            plotlayout_before = plotlayout
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
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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
            plotdata_before = plotdata
            plotlayout_before = plotlayout
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
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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
            
            if !haskey(registry, selected_folder_main)
                msg = "Dataset '$(selected_folder_main)' not found in registry."
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

            y_positive = yCoord < 0 ? abs(yCoord) : yCoord
            plotdata, plotlayout, xSpectraMz, ySpectraMz = xySpectrumPlot(msi_data, xCoord, y_positive, imgWidth, imgHeight, selected_folder_main, mask_path=mask_path_for_plot)
            plotdata_before = plotdata
            plotlayout_before = plotlayout
            
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
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
            end
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