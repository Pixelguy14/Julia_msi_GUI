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
using .MSI_src: MSIData, OpenMSIData, process_spectrum, IterateSpectra, ImzMLSource, _iterate_spectra_fast, MzMLSource, find_mass, ViridisPalette, get_mz_slice, get_multiple_mz_slices, quantize_intensity, save_bitmap, median_filter, save_bitmap, downsample_spectrum, TrIQ, precompute_analytics, ImportMzmlFile, generate_colorbar_image, load_and_prepare_mask, set_global_mz_range!, main_precalculation, MutableSpectrum, execute_full_preprocessing

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

function validate_parse(validation_errors::Vector{String}, param_str::String, param_name::String, target_type::Type, step_name::String)
    println("DEBUG: Validating ($step_name) Parameter '$param_name'. Received value: '$param_str'")
    if isempty(param_str)
        push!(validation_errors, "($step_name) Parameter '$param_name' is empty.")
        return nothing
    end
    val = tryparse(target_type, param_str)
    if val === nothing
        push!(validation_errors, "($step_name) Parameter '$param_name' ('$param_str') is not a valid $(target_type).")
        return nothing
    end
    return val
end

# Helper function to check if a pipeline step is enabled
function is_step_enabled(step_name::String, pipeline_order::Vector{Dict{String, Any}})
    for step in pipeline_order
        if get(step, "name", "") == step_name
            return get(step, "enabled", false)
        end
    end
    return false # Default to disabled if step not found
end

function get_processed_mean_spectrum(spectra::Vector{MutableSpectrum}; num_bins=2000)
    # 1. Find global m/z range from all spectra
    min_mz, max_mz = Inf, -Inf
    for s in spectra
        if !isempty(s.mz)
            min_mz = min(min_mz, minimum(s.mz))
            max_mz = max(max_mz, maximum(s.mz))
        end
    end

    if !isfinite(min_mz)
        return Float64[], Float64[]
    end

    # 2. Create bins
    mz_bins = range(min_mz, stop=max_mz, length=num_bins)
    intensity_sum = zeros(Float64, num_bins)
    bin_step = step(mz_bins)
    inv_bin_step = 1.0 / bin_step

    # 3. Bin intensities
    for s in spectra
        for i in eachindex(s.mz)
            bin_index = trunc(Int, (s.mz[i] - min_mz) * inv_bin_step + 1.0)
            final_index = clamp(bin_index, 1, num_bins)
            intensity_sum[final_index] += s.intensity[i]
        end
    end

    # 4. Average and return
    if isempty(spectra)
        return collect(mz_bins), intensity_sum
    end
    average_intensity = intensity_sum ./ length(spectra)
    return collect(mz_bins), average_intensity
end

function get_processed_sum_spectrum(spectra::Vector{MutableSpectrum}; num_bins=2000)
    min_mz, max_mz = Inf, -Inf
    for s in spectra
        if !isempty(s.mz)
            min_mz = min(min_mz, minimum(s.mz))
            max_mz = max(max_mz, maximum(s.mz))
        end
    end

    if !isfinite(min_mz)
        return Float64[], Float64[]
    end

    mz_bins = range(min_mz, stop=max_mz, length=num_bins)
    intensity_sum = zeros(Float64, num_bins)
    bin_step = step(mz_bins)
    inv_bin_step = 1.0 / bin_step

    for s in spectra
        for i in eachindex(s.mz)
            bin_index = trunc(Int, (s.mz[i] - min_mz) * inv_bin_step + 1.0)
            final_index = clamp(bin_index, 1, num_bins)
            intensity_sum[final_index] += s.intensity[i]
        end
    end
    
    return collect(mz_bins), intensity_sum
end

@genietools

# == Reactive code ==
# Reactive code to make the UI interactive
@app begin
    # == Loading Screen Variables ==
    @in is_initializing = true
    @in initialization_message = "Initializing..."
    # Loading animations and readonly / disable elements are all handled by this variable.
    @in is_processing = false

    # == SLICE GENERATOR TAB VARIABLES ==
    # File selection and batch processing
    @in file_route=""
    @in file_name=""
    @in btnSearch=false # To search for files in your device
    @in btnAddBatch = false
    @in clear_batch_btn = false
    @out batch_file_count = 0
    @in selected_files = String[]
    @out full_route="" # Saves the route where imzML and mzML files are located

    # Mass-to-charge parameters
    @in Nmass="0.0" # Mass-to-charge ratio(s) of interest
    @in Tol=0.1 # Mass-to-charge ratio tolerance
    @in colorLevel=20 # Color levels for visualization

    # Processing toggles
    @in triqEnabled=false # Threshold Intensity Quantization
    @in MFilterEnabled=false # Median Filter
    @in maskEnabled=false # Use Mask To Filter Data
    @in triqProb=0.98 # TrIQ probability parameter

    # Spectrum selection and coordinates
    @in idSpectrum=0 # Spectrum ID for ID-based plots
    @in xCoord=0 # X coordinate for coordinate-based plots
    @in yCoord=0 # Y coordinate for coordinate-based plots
    @in SpectraEnabled=false # Enables xCoord and yCoord inputs when spectral data is loaded

    # Plot generation triggers
    @in mainProcess=false # To generate images/slices
    @in createMeanPlot=false # To generate mean spectrum plot
    @in createXYPlot=false # To generate spectrum plot according to xy values
    @in createNSpectrumPlot=false # To generate spectrum plot according to spectrum order
    @in createSumPlot=false # To generate sum of all spectrum plots
    @in image3dPlot=false # To generate 3d plot based on current image
    @in triq3dPlot=false # To generate 3d plot based on current triq image
    @in imageCPlot=false # To generate contour plots of current image
    @in triqCPlot=false # To generate contour plots of current triq image

    # Image navigation controls
    @in imgPlus=false # Next image in normal mode
    @in imgMinus=false # Previous image in normal mode
    @in imgPlusT=false # Next image in TrIQ mode
    @in imgMinusT=false # Previous image in TrIQ mode

    # Optical image overlay
    @in imgTrans=1.0 # Transparency level for optical overlay
    @in btnOptical=false # Load optical image over normal image
    @in btnOpticalT=false # Load optical image over TrIQ image
    @in opticalOverTriq=false # Toggle optical overlay mode

    # Messages and status
    @out msg="" # Main status message
    @out msgimg="" # Image status message
    @out msgtriq="" # TrIQ status message

    # == CONVERTER TAB VARIABLES ==
    @in left_tab = "generator" # Active left tab (generator, converter, pre_treatment)
    @out mzml_full_route = "" # Path to .mzML file
    @out sync_full_route = "" # Path to .txt synchronization file
    @in btnSearchMzml = false # Trigger mzML file search
    @in btnSearchSync = false # Trigger sync file search
    @in convert_process = false # Start conversion process
    @out progress_conversion = false # Conversion progress indicator
    @out msg_conversion = "" # Conversion status message
    @out btnConvertDisable = true # Disable convert button when files not selected

    # == PRE-TREATMENT TAB VARIABLES ==
    # File selection and batch
    @in pre_tab = "stabilization" # Active preprocessing subtab

    # Subset processing
    @in enable_subset_processing = false # Enable processing only first N spectra
    @in spectra_subset_size = 100 # Number of spectra for subset processing

    # Internal standards management
    @in enable_standards = true # Use internal standards for calibration
    @in reference_peaks_list = [
        Dict("mz" => 137.0244, "label" => "DHB_fragment"),
        Dict("mz" => 155.0349, "label" => "DHB_M+H"),
    ]
    @in addReferencePeak = false # Add new reference peak
    @in remove_peak_trigger = false # Remove reference peak
    @in export_standards_btn = false # Export standards to JSON
    @in import_standards_btn = false # Import standards from JSON

    # Pipeline step management
    @in pipeline_step_order = [
        Dict("name" => "stabilization", "label" => "Stabilization", "enabled" => true),
        Dict("name" => "smoothing", "label" => "Smoothing", "enabled" => true),
        Dict("name" => "baseline_correction", "label" => "Baseline Correction", "enabled" => true),
        Dict("name" => "peak_picking", "label" => "Peak Picking", "enabled" => true),
        Dict("name" => "peak_selection", "label" => "Peak Selection", "enabled" => true),
        Dict("name" => "calibration", "label" => "Calibration", "enabled" => true),
        Dict("name" => "peak_alignment", "label" => "Peak Alignment", "enabled" => true),
        Dict("name" => "normalization", "label" => "Normalization", "enabled" => true),
        Dict("name" => "peak_binning", "label" => "Peak Binning", "enabled" => true)
    ]
    @in action_index = -1 # Index for step operations
    @in move_step_up_trigger = false # Move step up in pipeline
    @in move_step_down_trigger = false # Move step down in pipeline
    @in toggle_step_trigger = false # Toggle step enabled/disabled
    @out current_pipeline_step = "" # Current running step in full pipeline

    # Preprocessing method parameters
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

    # Suggested parameter values
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

    # Pipeline control triggers
    @in run_full_pipeline = false # Trigger full pipeline execution
    @in recalculate_suggestions_btn = false # Recalculate parameter suggestions
    @in export_params_btn = false # Export parameters to file
    @in import_params_btn = false # Import parameters from file
    @in save_feature_matrix_btn = false # Save feature matrix results

    # Preprocessing results
    @in selected_spectrum_id_for_plot = 1
    @in last_plot_type = "single"
    @in feature_matrix_result::Union{Nothing, Matrix{Float64}} = nothing
    @in bin_info_result::Union{Nothing, Vector} = nothing

    # == RIGHT PANEL VARIABLES (intDivStyle-right) ==
    # Tab management
    @out tabIDs=["tab0","tab1","tab2","tab3","tab4"]
    @out tabLabels=["Image", "TrIQ", "Spectrum Plot", "Topography Plot","Surface Plot"]
    @in selectedTab="tab0"

    # Compare dialog tabs
    @out CompTabIDsLeft=["tab0","tab1","tab2","tab3","tab4"]
    @out CompTabLabelsLeft=["Image", "TrIQ", "Spectrum Plot", "Topography Plot","Surface Plot"]
    @in CompSelectedTabLeft="tab0"
    @out CompTabIDsRight=["tab0","tab1","tab2","tab3","tab4"]
    @out CompTabLabelsRight=["Image", "TrIQ", "Spectrum Plot", "Topography Plot","Surface Plot"]
    @in CompSelectedTabRight="tab0"

    # Compare dialog controls
    @in CompareDialog=false
    @in compareBtn=false # Open compare dialog
    @in imgPlusCompLeft=false # Next image in compare left panel
    @in imgMinusCompLeft=false # Previous image in compare left panel
    @in imgPlusTCompLeft=false # Next TrIQ image in compare left panel
    @in imgMinusTCompLeft=false # Previous TrIQ image in compare left panel
    @in imgPlusCompRight=false # Next image in compare right panel
    @in imgMinusCompRight=false # Previous image in compare right panel
    @in imgPlusTCompRight=false # Next TrIQ image in compare right panel
    @in imgMinusTCompRight=false # Previous TrIQ image in compare right panel

    # Image display variables
    @out imgInt="/.bmp" # Normal image interface
    @out imgIntT="/.bmp" # TrIQ image interface
    @out colorbar="/.png" # Normal colorbar
    @out colorbarT="/.png" # TrIQ colorbar

    # Compare dialog images
    @out imgIntCompLeft="/.bmp" # Left compare normal image
    @out imgIntTCompLeft="/.bmp" # Left compare TrIQ image
    @out colorbarCompLeft="/.png" # Left compare normal colorbar
    @out colorbarTCompLeft="/.png" # Left compare TrIQ colorbar
    @out imgIntCompRight="/.bmp" # Right compare normal image
    @out imgIntTCompRight="/.bmp" # Right compare TrIQ image
    @out colorbarCompRight="/.png" # Right compare normal colorbar
    @out colorbarTCompRight="/.png" # Right compare TrIQ colorbar

    @out imgWidth=0
    @out imgHeight=0

    # Compare dialog messages
    @out msgimgCompLeft=""
    @out msgtriqCompLeft=""
    @out msgimgCompRight=""
    @out msgtriqCompRight=""

    # == BATCH PROCESSING & REGISTRY VARIABLES ==
    @private registry_init_done = false
    @in refetch_folders = false
    @in available_folders = String[]
    @in image_available_folders = String[]
    @out registry_path = abspath(joinpath(@__DIR__, "public", "registry.json"))

    # Folder selection state
    @in selected_folder_main = ""
    @in selected_folder_compare_left = ""
    @in selected_folder_compare_right = ""

    # Progress reporting
    @out overall_progress = 0.0
    @out progress_message = ""

    # Batch summary
    @in showBatchSummary = false
    @out batch_summary = ""

    # == METADATA VARIABLES ==
    @in showMetadataDialog = false
    @in showMetadataBtn = false
    @out metadata_columns = []
    @out metadata_rows = []
    @out btnMetadataDisable = false
    @in selected_folder_metadata = ""

    # == DATA MANAGEMENT VARIABLES ==
    # Centralized MSIData object
    @out msi_data::Union{MSIData, Nothing} = nothing

    # Image file management
    @out text_nmass="" # For specific mass charge image creation
    @in msi_bmp=sort(filter(filename -> startswith(filename, "MSI_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_msi_png=sort(filter(filename -> startswith(filename, "colorbar_MSI_") && endswith(filename, ".png"), readdir("public")),lt=natural)
    @in triq_bmp=sort(filter(filename -> startswith(filename, "TrIQ_") && endswith(filename, ".bmp"), readdir("public")),lt=natural)
    @in col_triq_png=sort(filter(filename -> startswith(filename, "colorbar_TrIQ_") && endswith(filename, ".png"), readdir("public")),lt=natural)

    # Current image display
    @out current_msi=""
    @out current_col_msi=""
    @out current_triq=""
    @out current_col_triq=""
    @out current_msiCompLeft=""
    @out current_col_msiCompLeft=""
    @out current_triqCompLeft=""
    @out current_col_triqCompLeft=""
    @out current_msiCompRight=""
    @out current_col_msiCompRight=""
    @out current_triqCompRight=""
    @out current_col_triqCompRight=""

    # Optical image
    @out imgRoute=""

    # == PLOTTING VARIABLES ==
    # Image plots
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
    @out plotdataImgCompLeft=[traceImg]
    @out plotlayoutImgCompLeft=layoutImg
    @out plotdataImgCompRight=[traceImg]
    @out plotlayoutImgCompRight=layoutImg

    # TrIQ image plots
    @out plotdataImgT=[traceImg]
    @out plotlayoutImgT=layoutImg
    @out plotdataImgTCompLeft=[traceImg]
    @out plotlayoutImgTCompLeft=layoutImg
    @out plotdataImgTCompRight=[traceImg]
    @out plotlayoutImgTCompRight=layoutImg

    # Spectrum plots
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
        margin=attr(l=0,r=0,t=120,b=0,pad=0),
        legend=attr(
            x=1.0, 
            y=1.0,
            xanchor="right", 
            yanchor="top"
        )
    )
    traceSpectra=PlotlyBase.scatter(x=Vector{Float64}(), y=Vector{Float64}(), mode="lines", marker=attr(size=1, color="blue", opacity=0.1))
    @out plotdata=[traceSpectra]
    @out plotlayout=layoutSpectra

    # Preprocessing spectrum plots
    @out plotdata_before = [traceSpectra]
    @out plotlayout_before = layoutSpectra
    @out plotdata_after = [traceSpectra]
    @out plotlayout_after = layoutSpectra

    # Spectrum data
    @out xSpectraMz = Vector{Float64}()
    @out ySpectraMz = Vector{Float64}()

    # Contour plots
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
    traceContour=PlotlyBase.contour(x=Vector{Float64}(), y=Vector{Float64}(), mode="lines")
    @out plotdataC=[traceContour]
    @out plotlayoutC=layoutContour

    # 3D surface plots
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
    @out plotdata3d=[trace3D]
    @out plotlayout3d=layout3D

    # Interactive plot reactions
    @in data_click=Dict{String,Any}()

    # == TIME MEASUREMENT VARIABLES ==
    @out sTime=time()
    @out fTime=time()
    @out eTime=time()

    # == DIALOGS AND MESSAGES ==
    @in warning_msg=false

    # == Reactive handlers ==
    # Reactive handlers watch a variable and execute a block of code when its value changes
    # The onbutton handler will set the variable to false after the block is executed

    # This handler correctly uses pick_file and loads the selected file
    # as the active dataset for the UI.
    @onbutton btnSearch begin
        is_processing = true
        picked_route = pick_file(; filterlist="imzML,imzml,mzML,mzml")
        if isempty(picked_route)
            is_processing = false
            return
        end

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
                SpectraEnabled = true
                selected_folder_main = dataset_name
                
                # Update folder lists in UI
                all_folders = sort(collect(keys(registry)), lt=natural)
                img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)
                available_folders = deepcopy(all_folders)
                image_available_folders = deepcopy(img_folders)

                msg = "Successfully loaded pre-processed dataset: $(dataset_name)"
                
            end


            # --- Full Load Path ---
            
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

            if isempty(get(existing_entry, "metadata", Dict()))
                msg = "Performing first-time metadata analysis for: $(basename(picked_route))..."
                precompute_analytics(loaded_data)
            end
            
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

            SpectraEnabled = true
            
        catch e
            msi_data = nothing
            msg = "Error loading active file: $e"
            warning_msg = true
            SpectraEnabled = false
            btnMetadataDisable = true
            @error "File loading failed" exception=(e, catch_backtrace())
        finally
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
            end
            is_processing = false
        end
    end

    @onbutton export_params_btn begin
        is_processing = true
        params_to_export = Dict(
            "pipeline_step_order" => pipeline_step_order,
            "enable_standards" => enable_standards, # Export global flag
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
        
        # 1. Generate the JSON string
        json_string = JSON.json(params_to_export)

        # 2. ESCAPING (Crucial for stability)
        # We must escape backslashes (for Windows paths) and single quotes
        # so they don't break the JavaScript string literal.
        safe_json = replace(json_string, "\\" => "\\\\") 
        safe_json = replace(safe_json, "'" => "\\'")

        # 3. Create the JavaScript payload
        # We inject 'safe_json' into the JS 'encodeURIComponent'
        js_script = """
            var element = document.createElement('a');
            element.setAttribute('href', 'data:text/json;charset=utf-8,' + encodeURIComponent('$safe_json'));
            element.setAttribute('download', 'preprocessing_params.json');
            element.style.display = 'none';
            document.body.appendChild(element);
            element.click();
            document.body.removeChild(element);
        """

        # 4. Execute on the client
        run(__model__, js_script)
        is_processing = false
        msg = "Parameters exported."
    end

    @onbutton import_params_btn begin
        is_processing = true
        picked_file = pick_file(filterlist="json")
        if isempty(picked_file)
            is_processing = false
            return
        end

        try
            json_string = read(picked_file, String)
            params = JSON.parse(json_string)

            # This is a list of all known reactive variables that can be imported.
            # This prevents arbitrary variable assignment.
            known_params = [
                "stabilization_method", "smoothing_method", "smoothing_window", "smoothing_order",
                "baseline_method", "baseline_iterations", "baseline_window", "normalization_method",
                "alignment_method", "alignment_span", "alignment_tolerance", "alignment_tolerance_unit",
                "alignment_max_shift_ppm", "alignment_min_matched_peaks", "peak_picking_method",
                "peak_picking_snr_threshold", "peak_picking_half_window", "peak_picking_min_peak_prominence",
                "peak_picking_merge_peaks_tolerance", "peak_picking_min_peak_width_ppm",
                "peak_picking_max_peak_width_ppm", "peak_picking_min_peak_shape_r2", "binning_method",
                "binning_tolerance", "binning_tolerance_unit", "binning_frequency_threshold",
                "binning_min_peak_per_bin", "binning_max_bin_width_ppm", "binning_intensity_weighted_centers",
                "binning_num_uniform_bins", "calibration_fit_order", "calibration_ppm_tolerance",
                "peak_selection_min_snr", "peak_selection_min_fwhm_ppm", "peak_selection_max_fwhm_ppm",
                "peak_selection_min_shape_r2", "peak_selection_frequency_threshold", "peak_selection_correlation_threshold"
            ]

            for (key, value) in params
                if key == "reference_peaks_list"
                    reference_peaks_list = value
                elseif key == "pipeline_step_order"
                    pipeline_step_order = value
                elseif key == "enable_standards"
                    enable_standards = value
                elseif key in known_params
                    # Use getfield and setproperty! to update reactive variables by name
                    if hasfield(typeof(@__MODULE__), Symbol(key))
                        getfield(@__MODULE__, Symbol(key))[] = value
                    end
                else
                    @warn "Unknown parameter '$key' found in JSON file. Skipping."
                end
            end
            msg = "Parameters imported successfully from $(basename(picked_file))."
        catch e
            msg = "Failed to import parameters: $e"
            warning_msg = true
            @error "Parameter import failed" exception=(e, catch_backtrace())
        end
        is_processing = false
    end

    @onbutton export_standards_btn begin
        is_processing = true
        json_string = JSON.json(reference_peaks_list)
        safe_json = replace(json_string, "\\" => "\\\\")
        safe_json = replace(safe_json, "'" => "\\'")
        js_script = """
            var element = document.createElement('a');
            element.setAttribute('href', 'data:text/json;charset=utf-8,' + encodeURIComponent('$safe_json'));
            element.setAttribute('download', 'internal_standards.json');
            element.style.display = 'none';
            document.body.appendChild(element);
            element.click();
            document.body.removeChild(element);
        """
        run(__model__, js_script)
        is_processing = false
        msg = "Internal standards exported."
    end

    @onbutton import_standards_btn begin
        is_processing = true
        picked_file = pick_file(filterlist="json")
        if isempty(picked_file)
            return
        end

        try
            json_string = read(picked_file, String)
            new_standards = JSON.parse(json_string)
            # Basic validation
            if new_standards isa Vector && all(p -> p isa Dict && haskey(p, "mz") && haskey(p, "label"), new_standards)
                reference_peaks_list = new_standards
                msg = "Internal standards imported successfully from $(basename(picked_file))."
            else
                msg = "Invalid format for internal standards file."
                warning_msg = true
            end
        catch e
            msg = "Failed to import internal standards: $e"
            warning_msg = true
            @error "Standards import failed" exception=(e, catch_backtrace())
        end
        is_processing = false
    end

    @onbutton run_full_pipeline begin
        is_processing = true
        current_pipeline_step = "Initializing..."
        println("DEBUG: run_full_pipeline started.")
        # println("DEBUG: Current pipeline_step_order configuration: $pipeline_step_order")
        
        try
            # --- 1. Initial Checks and Data Loading ---
            println("DEBUG: Performing initial checks and data loading...")
            if msi_data === nothing || isempty(selected_folder_main)
                msg = "No dataset loaded. Please load a file using 'Select an imzMl / mzML file'."
                warning_msg = true
                println("DEBUG: $msg")
                is_processing = false
                return
            end

            registry = load_registry(registry_path)
            entry = get(registry, selected_folder_main, nothing)
            if entry === nothing
                msg = "Selected dataset '$(selected_folder_main)' not found in registry. Please reload the file."
                warning_msg = true
                println("DEBUG: $msg")
                return
            end
            target_path = entry["source_path"]

            # Ensure msi_data is for the currently selected file
            if full_route != target_path
                println("DEBUG: Active file path changed. Reloading MSI data: $(basename(target_path))")
                if msi_data !== nothing; close(msi_data); end
                msg = "Reloading $(basename(target_path)) for analysis..."
                full_route = target_path
                msi_data = OpenMSIData(target_path)
            else
                println("DEBUG: Using already loaded MSI data for $(basename(target_path)).")
            end

            # Mask path retrieval from registry
            local mask_path_for_pipeline::Union{String, Nothing} = nothing
            if maskEnabled
                println("DEBUG: Masking is ENABLED.")
                if get(entry, "has_mask", false)
                    mask_path_candidate = get(entry, "mask_path", "")
                    if isfile(mask_path_candidate)
                        mask_path_for_pipeline = mask_path_candidate
                        println("DEBUG: Using mask for pipeline: $(mask_path_for_pipeline)")
                    else
                        msg = "Mask enabled but file not found: $(mask_path_candidate). Aborting pipeline."
                        warning_msg = true
                        @warn msg
                        println("DEBUG: $msg")
                        return
                    end
                else
                    msg = "Mask enabled but no valid mask entry found for: $(selected_folder_main). Aborting pipeline."
                    warning_msg = true
                    @warn msg
                    println("DEBUG: $msg")
                    return
                end
            else
                println("DEBUG: Masking is DISABLED. No mask will be applied.")
            end

            # Apply mask if enabled to get indices to process
            spectrum_indices_to_process = collect(1:length(msi_data.spectra_metadata))
            if mask_path_for_pipeline !== nothing
                current_pipeline_step = "Applying mask..."
                println("DEBUG: Applying mask matrix to filter spectra...")
                mask_matrix = load_and_prepare_mask(mask_path_for_pipeline, msi_data.image_dims)
                masked_indices_set = get_masked_spectrum_indices(msi_data, mask_matrix)
                spectrum_indices_to_process = collect(masked_indices_set)
                if isempty(spectrum_indices_to_process)
                    msg = "No spectra remaining after applying mask. Aborting pipeline."
                    warning_msg = true
                    println("DEBUG: $msg")
                    return
                end
                println("DEBUG: $(length(spectrum_indices_to_process)) spectra remaining after mask application.")
            else
                println("DEBUG: No mask applied. Processing all $(length(msi_data.spectra_metadata)) spectra.")
            end

            # Apply subset processing if enabled
            if enable_subset_processing && spectra_subset_size > 0
                n_total = length(spectrum_indices_to_process)
                n_to_process = min(spectra_subset_size, n_total)
                spectrum_indices_to_process = spectrum_indices_to_process[1:n_to_process]
                println("DEBUG: Subset processing enabled. Processing first $(length(spectrum_indices_to_process)) of $n_total spectra.")
            end

            # Initialize spectra data structure
            current_pipeline_step = "Loading spectra..."
            println("DEBUG: Loading $(length(spectrum_indices_to_process)) spectra into MutableSpectrum objects...")
            current_spectra = Vector{MutableSpectrum}(undef, length(spectrum_indices_to_process))
            
            Threads.@threads for i in 1:length(spectrum_indices_to_process)
                original_idx = spectrum_indices_to_process[i]
                mz, intensity = GetSpectrum(msi_data, original_idx) # Fetch mz and intensity for the current spectrum
                current_spectra[i] = MutableSpectrum(original_idx, copy(Float64.(mz)), copy(Float64.(intensity)), NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), NTuple{6, Float64}}[])
            end
            println("DEBUG: All spectra loaded into temporary structure for processing.")


            # --- 2. Parameter Assembly with Validation ---
            current_pipeline_step = "Configuring parameters..."
            println("DEBUG: Configuring parameters and validating enabled steps...")
            ref_peaks = Dict{Float64, String}(p["mz"] => p["label"] for p in reference_peaks_list)

            final_params = Dict{Symbol, Dict{Symbol, Any}}()
            validation_errors = String[]
            
            # --- Stabilization ---
            println("DEBUG: Checking Stabilization step (name: stabilization)")
            if is_step_enabled("stabilization", pipeline_step_order)
                println("DEBUG: Stabilization step is ENABLED. Setting method: $(stabilization_method).")
                final_params[:Stabilization] = Dict{Symbol, Any}(:method => Symbol(stabilization_method))
            else
                println("DEBUG: Stabilization step is DISABLED. Skipping.")
            end

            # --- Smoothing ---
            println("DEBUG: Checking Smoothing step (name: smoothing)")
            if is_step_enabled("smoothing", pipeline_step_order)
                println("DEBUG: Smoothing step is ENABLED. Validating parameters.")
                window_val = validate_parse(validation_errors, smoothing_window, "Window", Int, "Smoothing")
                order_val = validate_parse(validation_errors, smoothing_order, "Order", Int, "Smoothing")
                final_params[:Smoothing] = Dict{Symbol, Any}(
                    :method => Symbol(smoothing_method),
                    :window => something(window_val, 9),
                    :order => something(order_val, 2)
                )
                if window_val !== nothing && window_val < 1
                    push!(validation_errors, "(Smoothing) Window must be positive.")
                end
                if order_val !== nothing && order_val < 0
                    push!(validation_errors, "(Smoothing) Order must be non-negative.")
                end
                println("DEBUG: Smoothing parameters set: method=$(smoothing_method), window=$(something(window_val, 9)), order=$(something(order_val, 2)).")
            else
                println("DEBUG: Smoothing step is DISABLED. Skipping parameter validation.")
            end

            # --- Baseline Correction ---
            println("DEBUG: Checking Baseline Correction step (name: baseline_correction)")
            if is_step_enabled("baseline_correction", pipeline_step_order)
                println("DEBUG: Baseline Correction step is ENABLED. Validating parameters.")
                iterations_val = validate_parse(validation_errors, baseline_iterations, "Iterations", Int, "Baseline Correction")
                baseline_window_val = validate_parse(validation_errors, baseline_window, "Window", Int, "Baseline Correction")
                final_params[:BaselineCorrection] = Dict{Symbol, Any}(
                    :method => Symbol(baseline_method),
                    :iterations => something(iterations_val, 100),
                    :window => something(baseline_window_val, 20)
                )
                if iterations_val !== nothing && iterations_val < 0
                    push!(validation_errors, "(Baseline Correction) Iterations must be non-negative.")
                end
                if baseline_window_val !== nothing && baseline_window_val < 1
                    push!(validation_errors, "(Baseline Correction) Window must be positive.")
                end
                println("DEBUG: Baseline Correction parameters set: method=$(baseline_method), iterations=$(something(iterations_val, 100)), window=$(something(baseline_window_val, 20)).")
            else
                println("DEBUG: Baseline Correction step is DISABLED. Skipping parameter validation.")
            end

            # --- Normalization ---
            println("DEBUG: Checking Normalization step (name: normalization)")
            if is_step_enabled("normalization", pipeline_step_order)
                println("DEBUG: Normalization step is ENABLED. Setting method: $(normalization_method).")
                final_params[:Normalization] = Dict{Symbol, Any}(:method => Symbol(normalization_method))
            else
                println("DEBUG: Normalization step is DISABLED. Skipping.")
            end


            # --- Peak Picking ---
            println("DEBUG: Checking Peak Picking step (name: peak_picking)")
            if is_step_enabled("peak_picking", pipeline_step_order)
                println("DEBUG: Peak Picking step is ENABLED. Validating parameters.")
                snr_threshold_val = validate_parse(validation_errors, peak_picking_snr_threshold, "SNR Threshold", Float64, "Peak Picking")
                half_window_val = validate_parse(validation_errors, peak_picking_half_window, "Half Window", Int, "Peak Picking")
                min_peak_prominence_val = validate_parse(validation_errors, peak_picking_min_peak_prominence, "Min Prominence", Float64, "Peak Picking")
                merge_peaks_tolerance_val = validate_parse(validation_errors, peak_picking_merge_peaks_tolerance, "Merge Tolerance", Float64, "Peak Picking")
                final_params[:PeakPicking] = Dict{Symbol, Any}(
                    :method => Symbol(peak_picking_method),
                    :snr_threshold => something(snr_threshold_val, 3.0),
                    :half_window => something(half_window_val, 10),
                    :min_peak_prominence => something(min_peak_prominence_val, 0.1),
                    :merge_peaks_tolerance => something(merge_peaks_tolerance_val, 0.002)
                )
                if snr_threshold_val !== nothing && snr_threshold_val < 0
                    push!(validation_errors, "(Peak Picking) SNR Threshold must be non-negative.")
                end
                if half_window_val !== nothing && half_window_val < 1
                    push!(validation_errors, "(Peak Picking) Half Window must be positive.")
                end
                if min_peak_prominence_val !== nothing && (min_peak_prominence_val < 0 || min_peak_prominence_val > 1)
                    push!(validation_errors, "(Peak Picking) Min Prominence must be between 0 and 1.")
                end
                if merge_peaks_tolerance_val !== nothing && merge_peaks_tolerance_val < 0
                    push!(validation_errors, "(Peak Picking) Merge Tolerance must be non-negative.")
                end
                println("DEBUG: Peak Picking parameters set: method=$(peak_picking_method), snr_threshold=$(something(snr_threshold_val, 3.0)), half_window=$(something(half_window_val, 10))...")
            else
                println("DEBUG: Peak Picking step is DISABLED. Skipping parameter validation.")
            end

            # --- Peak Selection ---
            println("DEBUG: Checking Peak Selection step (name: peak_selection)")
            if is_step_enabled("peak_selection", pipeline_step_order)
                println("DEBUG: Peak Selection step is ENABLED. Validating parameters.")
                min_snr_val = validate_parse(validation_errors, peak_selection_min_snr, "Min SNR", Float64, "Peak Selection")
                min_fwhm_ppm_val = validate_parse(validation_errors, peak_selection_min_fwhm_ppm, "Min FWHM", Float64, "Peak Selection")
                max_fwhm_ppm_val = validate_parse(validation_errors, peak_selection_max_fwhm_ppm, "Max FWHM", Float64, "Peak Selection")
                min_shape_r2_val = validate_parse(validation_errors, peak_selection_min_shape_r2, "Min Shape R2", Float64, "Peak Selection")
                final_params[:PeakSelection] = Dict{Symbol, Any}(
                    :min_snr => something(min_snr_val, 0.0),
                    :min_fwhm_ppm => something(min_fwhm_ppm_val, 0.0),
                    :max_fwhm_ppm => something(max_fwhm_ppm_val, Inf),
                    :min_shape_r2 => something(min_shape_r2_val, 0.0)
                )
                if min_snr_val !== nothing && min_snr_val < 0
                    push!(validation_errors, "(Peak Selection) Min SNR must be non-negative.")
                end
                if min_fwhm_ppm_val !== nothing && min_fwhm_ppm_val < 0
                    push!(validation_errors, "(Peak Selection) Min FWHM must be non-negative.")
                end
                if max_fwhm_ppm_val !== nothing && max_fwhm_ppm_val < 0
                    push!(validation_errors, "(Peak Selection) Max FWHM must be non-negative.")
                end
                if min_shape_r2_val !== nothing && (min_shape_r2_val < 0 || min_shape_r2_val > 1)
                    push!(validation_errors, "(Peak Selection) Min Shape R2 must be between 0 and 1.")
                end
                println("DEBUG: Peak Selection parameters set: min_snr=$(something(min_snr_val, 0.0)), min_fwhm_ppm=$(something(min_fwhm_ppm_val, 0.0))...")
            else
                println("DEBUG: Peak Selection step is DISABLED. Skipping parameter validation.")
            end


            # --- Calibration ---
            println("DEBUG: Checking Calibration step (name: calibration)")
            if is_step_enabled("calibration", pipeline_step_order)
                println("DEBUG: Calibration step is ENABLED. Validating parameters.")
                ppm_tolerance_cal_val = validate_parse(validation_errors, calibration_ppm_tolerance, "PPM Tolerance", Float64, "Calibration")
                fit_order_val = validate_parse(validation_errors, calibration_fit_order, "Fit Order", Int, "Calibration")
                final_params[:Calibration] = Dict{Symbol, Any}(
                    :method => :internal_standards, # Fixed method
                    :ppm_tolerance => something(ppm_tolerance_cal_val, 20.0),
                    :fit_order => something(fit_order_val, 1) # Default to linear
                )
                if ppm_tolerance_cal_val !== nothing && ppm_tolerance_cal_val < 0
                    push!(validation_errors, "(Calibration) PPM Tolerance must be non-negative.")
                end
                if fit_order_val !== nothing && (fit_order_val < 0 || fit_order_val > 2)
                    push!(validation_errors, "(Calibration) Fit Order must be 0, 1, or 2.")
                end
                if enable_standards && isempty(ref_peaks)
                    push!(validation_errors, "(Calibration) Internal Standards are enabled, but no reference peaks are defined.")
                end
                println("DEBUG: Calibration parameters set: ppm_tolerance=$(something(ppm_tolerance_cal_val, 20.0)), fit_order=$(something(fit_order_val, 1)).")
            else
                println("DEBUG: Calibration step is DISABLED. Skipping parameter validation.")
            end


            # --- Peak Alignment ---
            println("DEBUG: Checking Peak Alignment step (name: peak_alignment)")
            if is_step_enabled("peak_alignment", pipeline_step_order)
                println("DEBUG: Peak Alignment step is ENABLED. Validating parameters.")
                alignment_tolerance_val = validate_parse(validation_errors, alignment_tolerance, "Tolerance", Float64, "Peak Alignment")
                final_params[:PeakAlignment] = Dict{Symbol, Any}(
                    :method => Symbol(alignment_method),
                    :tolerance => something(alignment_tolerance_val, 0.002),
                    :tolerance_unit => Symbol(alignment_tolerance_unit)
                )
                if alignment_tolerance_val !== nothing && alignment_tolerance_val < 0
                    push!(validation_errors, "(Peak Alignment) Tolerance must be non-negative.")
                end
                println("DEBUG: Peak Alignment parameters set: method=$(alignment_method), tolerance=$(something(alignment_tolerance_val, 0.002)), tolerance_unit=$(alignment_tolerance_unit).")
            else
                println("DEBUG: Peak Alignment step is DISABLED. Skipping parameter validation.")
            end


            # --- Peak Binning ---
            println("DEBUG: Checking Peak Binning step (name: peak_binning)")
            if is_step_enabled("peak_binning", pipeline_step_order)
                println("DEBUG: Peak Binning step is ENABLED. Validating parameters.")
                binning_tolerance_val = validate_parse(validation_errors, binning_tolerance, "Tolerance", Float64, "Peak Binning")
                min_peak_per_bin_val = validate_parse(validation_errors, binning_min_peak_per_bin, "Min Peaks Per Bin", Int, "Peak Binning")
                final_params[:PeakBinning] = Dict{Symbol, Any}(
                    :method => Symbol(binning_method),
                    :tolerance => something(binning_tolerance_val, 20.0),
                    :tolerance_unit => Symbol(binning_tolerance_unit),
                    :min_peak_per_bin => something(min_peak_per_bin_val, 3)
                )
                if binning_tolerance_val !== nothing && binning_tolerance_val < 0
                    push!(validation_errors, "(Peak Binning) Tolerance must be non-negative.")
                end
                if min_peak_per_bin_val !== nothing && min_peak_per_bin_val < 1
                    push!(validation_errors, "(Peak Binning) Min Peaks Per Bin must be positive.")
                end
                println("DEBUG: Peak Binning parameters set: method=$(binning_method), tolerance=$(something(binning_tolerance_val, 20.0)), min_peak_per_bin=$(something(min_peak_per_bin_val, 3))...")
            else
                println("DEBUG: Peak Binning step is DISABLED. Skipping parameter validation.")
            end


            if !isempty(validation_errors)
                msg = "Pipeline setup errors:\n" * join(validation_errors, "\n")
                warning_msg = true
                println("DEBUG: Validation errors encountered: $validation_errors")
                return
            end

            # Build pipeline steps from enabled steps in order
            pipeline_stp = [step["name"] for step in pipeline_step_order if step["enabled"]]
            println("DEBUG: Final enabled pipeline steps to execute: $pipeline_stp")
            
            # 3. Execute Pipeline
            current_pipeline_step = "Running preprocessing pipeline..."
            println("DEBUG: Starting pipeline execution with $(length(pipeline_stp)) enabled steps.")
            feature_matrix_result, bin_info_result = execute_full_preprocessing(
                current_spectra,
                final_params,
                pipeline_stp,
                ref_peaks,
                mask_path_for_pipeline
            ) do step
                current_pipeline_step = "Processing: $step"
                println("DEBUG: Processing step: $step")
            end
            println("DEBUG: Pipeline execution finished.")

            # 4. Update Results Display
            current_pipeline_step = "Updating results..."
            subset_label = enable_subset_processing ? " (from subset of $(length(current_spectra)) spectra)" : ""
            println("DEBUG: Updating results display after pipeline completion for plot type: $(last_plot_type)")

            if last_plot_type == "single"
                display_spectrum_idx = findfirst(s -> s.id == selected_spectrum_id_for_plot, current_spectra)
                if display_spectrum_idx !== nothing
                    processed_spectrum = current_spectra[display_spectrum_idx]
                    println("DEBUG: Displaying spectrum $(selected_spectrum_id_for_plot) after processing.")
                    
                    after_trace = PlotlyBase.scatter(
                        x=processed_spectrum.mz, 
                        y=processed_spectrum.intensity, 
                        mode="lines", 
                        name="Processed Spectrum"
                    )
                    
                    traces_after = [after_trace]

                    if !isempty(processed_spectrum.peaks)
                        peak_mzs = [p.mz for p in processed_spectrum.peaks]
                        peak_intensities = [p.intensity for p in processed_spectrum.peaks]
                        peak_trace = PlotlyBase.scatter(
                            x=peak_mzs, 
                            y=peak_intensities, 
                            mode="markers", 
                            name="Picked Peaks", 
                            marker=attr(color="red", size=8)
                        )
                        push!(traces_after, peak_trace)
                    end
                    
                    plotdata_after = traces_after
                    plotlayout_after = PlotlyBase.Layout(
                        title="After Preprocessing (Spectrum $(selected_spectrum_id_for_plot))$(subset_label)",
                        xaxis_title="m/z",
                        yaxis_title="Intensity",
                        legend=attr(x=0.98, y=0.98, xanchor="right", yanchor="top")
                    )
                else
                    println("DEBUG: Selected spectrum for display ($(selected_spectrum_id_for_plot)) not found in processed spectra.")
                end
            elseif last_plot_type == "mean"
                mz, intensity = get_processed_mean_spectrum(current_spectra)
                plotdata_after = [PlotlyBase.scatter(x=mz, y=intensity, mode="lines", name="Processed Mean Spectrum")]
                plotlayout_after = PlotlyBase.Layout(
                    title="After Preprocessing (Mean Spectrum)$(subset_label)",
                    xaxis_title="m/z",
                    yaxis_title="Average Intensity",
                    legend=attr(x=0.98, y=0.98, xanchor="right", yanchor="top")
                )
            elseif last_plot_type == "sum"
                mz, intensity = get_processed_sum_spectrum(current_spectra)
                plotdata_after = [PlotlyBase.scatter(x=mz, y=intensity, mode="lines", name="Processed Sum Spectrum")]
                plotlayout_after = PlotlyBase.Layout(
                    title="After Preprocessing (Sum Spectrum)$(subset_label)",
                    xaxis_title="m/z",
                    yaxis_title="Total Intensity",
                    legend=attr(x=0.98, y=0.98, xanchor="right", yanchor="top")
                )
            end

            # Save feature matrix if binning was performed
            if feature_matrix_result !== nothing
                output_dir = joinpath("public", selected_folder_main, "preprocessing_results")
                mkpath(output_dir)
                save_feature_matrix(feature_matrix_result, bin_info_result, output_dir)
                msg = "Pipeline completed successfully. Feature matrix saved."
                println("DEBUG: Feature matrix saved to $output_dir")
            else
                msg = "Pipeline completed successfully. No feature matrix generated (binning step not enabled)."
                println("DEBUG: $msg")
            end

        catch e
            msg = "Error during pipeline execution: $e"
            warning_msg = true
            @error "Pipeline failed" exception=(e, catch_backtrace())
            println("DEBUG: Pipeline caught an exception: $e")
        finally
            is_processing = false
            current_pipeline_step = ""
            println("DEBUG: run_full_pipeline finished (finally block).")
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
            end
        end
    end

    @onbutton recalculate_suggestions_btn begin
        is_processing = true
        if msi_data === nothing
            msg = "Please load a file first."
            warning_msg = true
            return
        end
        
        try
            msg = "Recalculating suggestions..."
            
            ref_peaks = Dict(p["mz"] => p["label"] for p in reference_peaks_list)
            recommended_params = main_precalculation(msi_data, reference_peaks=ref_peaks)
            

            for (step_name, params) in recommended_params
                for (param_key, value) in params
                    # Convert value to appropriate type before assignment
                    processed_value = if value === nothing
                        ""
                    elseif value isa Tuple
                        @warn "Skipping invalid parameter suggestion (tuple): $value for $param_key"
                        "" # Set to empty string for safety
                    elseif value isa Number
                        string(value)
                    else
                        string(value)
                    end
                    
                    if isempty(processed_value) && !(processed_value isa Number)
                        continue # Skip if processed_value is an empty string and not a number type
                    end

                    # Map recommended parameters to suggested_* reactive variables
                    if step_name == :Smoothing
                        if param_key == :window
                            suggested_smoothing_window = processed_value
                            smoothing_window = processed_value
                        elseif param_key == :order
                            suggested_smoothing_order = processed_value
                            smoothing_order = processed_value
                        end
                    elseif step_name == :BaselineCorrection
                        if param_key == :iterations
                            suggested_baseline_iterations = processed_value
                            baseline_iterations = processed_value
                        elseif param_key == :window
                            suggested_baseline_window = processed_value
                            baseline_window = processed_value
                        end
                    elseif step_name == :PeakAlignment
                        if param_key == :span
                            suggested_alignment_span = processed_value
                            alignment_span = processed_value
                        elseif param_key == :tolerance
                            suggested_alignment_tolerance = processed_value
                            alignment_tolerance = processed_value
                        elseif param_key == :max_shift_ppm
                            suggested_alignment_max_shift_ppm = processed_value
                            alignment_max_shift_ppm = processed_value
                        elseif param_key == :min_matched_peaks
                            suggested_alignment_min_matched_peaks = processed_value
                            alignment_min_matched_peaks = processed_value
                        end
                    elseif step_name == :Calibration
                        if param_key == :fit_order
                            suggested_calibration_fit_order = processed_value
                            calibration_fit_order = processed_value
                        elseif param_key == :ppm_tolerance
                            suggested_calibration_ppm_tolerance = processed_value
                            calibration_ppm_tolerance = processed_value
                        end
                    elseif step_name == :PeakPicking
                        if param_key == :snr_threshold
                            suggested_peak_picking_snr_threshold = processed_value
                            peak_picking_snr_threshold = processed_value
                        elseif param_key == :half_window
                            suggested_peak_picking_half_window = processed_value
                            peak_picking_half_window = processed_value
                        elseif param_key == :min_peak_prominence
                            suggested_peak_picking_min_peak_prominence = processed_value
                            peak_picking_min_peak_prominence = processed_value
                        elseif param_key == :merge_peaks_tolerance
                            suggested_peak_picking_merge_peaks_tolerance = processed_value
                            peak_picking_merge_peaks_tolerance = processed_value
                        elseif param_key == :min_peak_width_ppm
                            suggested_peak_picking_min_peak_width_ppm = processed_value
                            peak_picking_min_peak_width_ppm = processed_value
                        elseif param_key == :max_peak_width_ppm
                            suggested_peak_picking_max_peak_width_ppm = processed_value
                            peak_picking_max_peak_width_ppm = processed_value
                        elseif param_key == :min_peak_shape_r2
                            suggested_peak_picking_min_peak_shape_r2 = processed_value
                            peak_picking_min_peak_shape_r2 = processed_value
                        end
                    elseif step_name == :PeakSelection
                        if param_key == :min_snr
                            suggested_peak_selection_min_snr = processed_value
                            peak_selection_min_snr = processed_value
                        elseif param_key == :min_fwhm_ppm
                            suggested_peak_selection_min_fwhm_ppm = processed_value
                            peak_selection_min_fwhm_ppm = processed_value
                        elseif param_key == :max_fwhm_ppm
                            suggested_peak_selection_max_fwhm_ppm = processed_value
                            peak_selection_max_fwhm_ppm = processed_value
                        elseif param_key == :min_shape_r2
                            suggested_peak_selection_min_shape_r2 = processed_value
                            peak_selection_min_shape_r2 = processed_value
                        elseif param_key == :frequency_threshold
                            suggested_peak_selection_frequency_threshold = processed_value
                            peak_selection_frequency_threshold = processed_value
                        elseif param_key == :correlation_threshold
                            suggested_peak_selection_correlation_threshold = processed_value
                            peak_selection_correlation_threshold = processed_value
                        end
                    elseif step_name == :PeakBinning
                        if param_key == :tolerance
                            suggested_binning_tolerance = processed_value
                            binning_tolerance = processed_value
                        elseif param_key == :frequency_threshold
                            suggested_binning_frequency_threshold = processed_value
                            binning_frequency_threshold = processed_value
                        elseif param_key == :min_peak_per_bin
                            suggested_binning_min_peak_per_bin = processed_value
                            binning_min_peak_per_bin = processed_value
                        elseif param_key == :max_bin_width_ppm
                            suggested_binning_max_bin_width_ppm = processed_value
                            binning_max_bin_width_ppm = processed_value
                        elseif param_key == :num_uniform_bins
                            suggested_binning_num_uniform_bins = processed_value
                            binning_num_uniform_bins = processed_value
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
            
            msg = "Suggestions have been recalculated."
        catch e
            msg = "Failed to recalculate suggestions: $e"
            warning_msg = true
            @error "Recalculation failed" exception=(e, catch_backtrace())
        finally
            is_processing = false
        end
    end

    @onbutton addReferencePeak begin
        is_processing = true
        new_list = deepcopy(reference_peaks_list)
        push!(new_list, Dict("mz" => 0.0, "label" => ""))
        reference_peaks_list = new_list # Assign new list to trigger reactivity
        is_processing = false
    end

    @onbutton remove_peak_trigger begin
        is_processing = true
        if action_index > -1
            julia_index = action_index + 1
            new_list = deepcopy(reference_peaks_list)
            if 1 <= julia_index <= length(new_list)
                deleteat!(new_list, julia_index)
                reference_peaks_list = new_list
            end
            action_index = -1 # Reset
        end
        is_processing = false
    end

    @onbutton move_step_up_trigger begin
        is_processing = true
        if action_index > -1
            julia_index = action_index + 1
            if julia_index > 1
                new_order = deepcopy(pipeline_step_order)
                temp = new_order[julia_index]
                new_order[julia_index] = new_order[julia_index - 1]
                new_order[julia_index - 1] = temp
                pipeline_step_order = new_order
            end
            action_index = -1 # Reset
        end
        is_processing = false
    end

    @onbutton move_step_down_trigger begin
        is_processing = true
        if action_index > -1
            julia_index = action_index + 1
            if julia_index < length(pipeline_step_order)
                new_order = deepcopy(pipeline_step_order)
                temp = new_order[julia_index]
                new_order[julia_index] = new_order[julia_index + 1]
                new_order[julia_index + 1] = temp
                pipeline_step_order = new_order
            end
            action_index = -1 # Reset
        end
        is_processing = false
    end
    
    @onbutton toggle_step_trigger begin
        is_processing = true
        if action_index > -1
            julia_index = action_index + 1
            if 1 <= julia_index <= length(pipeline_step_order)
                new_order = deepcopy(pipeline_step_order)
                new_order[julia_index]["enabled"] = !new_order[julia_index]["enabled"]
                pipeline_step_order = new_order
            end
            action_index = -1 # Reset
        end
        is_processing = false
    end

    # This new handler correctly adds the file from full_route to the batch list.
    @onbutton btnAddBatch begin
        is_processing = true
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
        is_processing = false
    end

    @onbutton clear_batch_btn begin
        is_processing = true
        selected_files = String[]
        batch_file_count = 0
        msg = "Batch cleared"
        is_processing = false
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
        is_processing = true
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
        is_processing = false
    end

    @onbutton convert_process begin
        is_processing = true
        if isempty(mzml_full_route) || isempty(sync_full_route)
            msg_conversion = "Please select both an .mzML file and a .txt sync file."
            warning_msg = true
            return
        end

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
            is_processing = false
            # Re-enable button if files are still selected
            btnConvertDisable = isempty(mzml_full_route) || isempty(sync_full_route)
        end
    end
    
    @onbutton mainProcess @time begin
        # --- UI State Update ---
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
            is_processing = true

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
            is_processing = false
            SpectraEnabled = true
            overall_progress = 0.0
            #println("Done")
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

        is_processing = true

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
            last_plot_type = "mean"
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
                    is_processing = false
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

        is_processing = true
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
            last_plot_type = "sum"
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
            is_processing = false
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

        is_processing = true
        msg = "Loading plot for $(selected_folder_main)..."

        try
            sTime = time()
            registry = load_registry(registry_path)
            
            # Add error handling for registry access
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

            # Convert to positive coordinates for processing
            y_positive = yCoord < 0 ? abs(yCoord) : yCoord
            plotdata, plotlayout, xSpectraMz, ySpectraMz, spectrum_id = xySpectrumPlot(msi_data, xCoord, y_positive, imgWidth, imgHeight, selected_folder_main, mask_path=mask_path_for_plot)
            plotdata_before = plotdata
            plotlayout_before = plotlayout
            last_plot_type = "single"
            selected_spectrum_id_for_plot = spectrum_id
            idSpectrum = spectrum_id # we set the same obtained spectrum id to the UI
            
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
            is_processing = false
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
            end
        end
    end

    @onbutton createNSpectrumPlot @time begin
        if isempty(selected_folder_main)
            msg = "No dataset selected. Please process a file and select a folder first."
            warning_msg = true
            return
        end

        is_processing = true
        msg = "Loading plot for $(selected_folder_main)..."

        try
            sTime = time()
            registry = load_registry(registry_path)
            
            # Add error handling for registry access
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

            # Call the new nSpectrumPlot function
            plotdata, plotlayout, xSpectraMz, ySpectraMz, spectrum_id = nSpectrumPlot(msi_data, idSpectrum, selected_folder_main, mask_path=mask_path_for_plot)
            plotdata_before = plotdata
            plotlayout_before = plotlayout
            last_plot_type = "single"
            selected_spectrum_id_for_plot = spectrum_id
            
            selectedTab = "tab2"
            fTime = time()
            eTime = round(fTime - sTime, digits=3)
            msg = "Plot loaded in $(eTime) seconds"
            log_memory_usage("nSpectrum Plot Generated", msi_data)
        catch e
            msg = "Could not retrieve spectrum: $e"
            warning_msg = true
            @error "nSpectrum plotting failed" exception=(e, catch_backtrace())
        finally
            is_processing = false
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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
        GC.gc() # Trigger garbage collection
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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
        is_processing = true

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
            is_processing = true
            SpectraEnabled=true
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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

        is_processing = true

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
            is_processing = true
            SpectraEnabled=true
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure Julia returns the freed memory to OS
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

        is_processing = true
        
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
            is_processing = true
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

        is_processing = true
        
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
            is_processing = true
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
            df = msi_data.spectrum_stats_df
            profile_count = 0
            if df !== nothing && hasproperty(df, :Mode)
                profile_count = count(==(MSI_src.PROFILE), df.Mode)
            end

            if profile_count > 0
                # Main spectrum trace
                traceSpectra = PlotlyBase.scatter(
                    x=xSpectraMz, 
                    y=ySpectraMz, 
                    marker=attr(size=1, color="blue", opacity=0.5), 
                    name="Spectrum", 
                    hoverinfo="x", 
                    hovertemplate="<b>m/z</b>: %{x:.4f}<extra></extra>", 
                    showlegend=false
                )
            else
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
            end
            
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
        is_processing = true
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
        is_processing = false
    end

    @onbutton btnOpticalT begin
        is_processing = true
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
        is_processing = false
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

    @onbutton refetch_folders begin
        is_processing = true
        # Re-load registry and update folder lists
        registry = load_registry(registry_path)
        all_folders = sort(collect(keys(registry)), lt=natural)
        img_folders = filter(folder -> get(get(registry, folder, Dict()), "is_imzML", false), all_folders)
        
        available_folders = deepcopy(all_folders)
        image_available_folders = deepcopy(img_folders)

        # For q-selects using image_available_folders
        if !isempty(image_available_folders)
            first_img_folder = first(image_available_folders)
            if isempty(selected_folder_main)
                selected_folder_main = first_img_folder
            end
            if isempty(selected_folder_compare_left)
                selected_folder_compare_left = first_img_folder
            end
            if isempty(selected_folder_compare_right)
                selected_folder_compare_right = first_img_folder
            end
        end

        # For q-selects using available_folders
        if !isempty(available_folders)
            if isempty(selected_folder_metadata)
                selected_folder_metadata = first(available_folders)
            end
        end
        is_processing = false
    end

    @mounted watchplots()

    @onchange isready @time begin
        # is_processing = true
        if isready && !registry_init_done
            initialization_message = "Pre-compiling functions at startup..."
            warmup_init()
            initialization_message = "Pre-compilation finished."
            try
                initialization_message = "Synchronizing registry with filesystem on backend init..."
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
                    initialization_message = "Registry changed, saving..."
                    save_registry(reg_path, registry)
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
                is_initializing = false # Hide loading screen when initialization is complete
            end
        end
        is_initializing = false # Hide loading screen when initialization is complete (current code is hidden due to incompatibility)
        initialization_message = "Done."
        log_memory_usage("App Ready", msi_data)
    end
    # is_processing = false
    GC.gc() # Trigger garbage collection
    if Sys.islinux()
        ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
    end
end
# == Pages ==
# Register a new route and the page that will be loaded on access
@page("/", "app.jl.html")
end