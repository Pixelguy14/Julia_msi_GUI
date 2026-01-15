# run_preprocessing.jl

using Printf
import Pkg
using CairoMakie
using DataFrames    # For creating dataframes
using CSV
using Statistics
using Interpolations

# --- Load the MSI_src Module ---
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src


# ===================================================================
# CONFIG: PLEASE FILL IN YOUR FILE PATHS HERE
# ===================================================================

# const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Thricoderma_etc/Imaging_interaccion_trichoderma_vs_streptomyces.imzML"
# const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/set de datos MS/Atropina_tuneo_fraq_20ev.mzML"
const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML"
# const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML"

const MASK_ROUTE = "/home/pixel/Documents/Cinvestav_2025/JuliaMSI/public/css/masks/Stomach_DHB_uncompressed.png"
# const MASK_ROUTE = "" 

const OUTPUT_DIR = "./test/results/preprocessing_results"

# Reference peaks for calibration and alignment
const reference_peaks = Dict(
    137.0244 => "DHB_fragment",
    155.0349 => "DHB_M+H", 
    177.0168 => "DHB_M+Na",
    496.3398 => "PC_16:0_16:0",
    520.3398 => "PC_16:0_18:1", 
    760.5851 => "PC_16:0_18:1_Na",
    391.2843 => "PDMS",
    413.2662 => "PDMS_Na",
    842.5092 => "Protein_standard",
    1045.532 => "Protein_standard",
    290.1747 => "Atropine [M+H]+", 
    304.1903 => "Scopolamine [M+H]+",
    124.0393 => "Tropine [M+H]+",
)

const PIPELINE_STP = [
    "stabilization",
    #"baseline_correction",
    "smoothing", 
    "peak_picking",
    "peak_selection",
    #"calibration",
    "peak_alignment",
    "normalization",
    "peak_binning"
]

# ===================================================================
# USER OVERRIDES: Manually specify parameters here
# ===================================================================
# This dictionary allows you to override any auto-detected parameters.
# The structure should match the output of `main_precalculation`.
# Example: Force a less aggressive peak prominence threshold.

const USER_OVERRIDES = Dict(
    :Stabilization => Dict(
        :method => :sqrt # Default stabilization method
    ),
    :PeakPicking => Dict(
        #:snr_threshold => 1.5,
        #:half_window => 5,
        :snr_threshold => 8.0,
        #:merge_peaks_tolerance => 2.5,
        #:half_window => 2
        #:half_window => 3
    ),
    # :Smoothing => Dict(
    #     :window => 11
    # )
    #:PeakAlignment => Dict(
    #    :tolerance => 0.01
    #)
    #:PeakSelection => Dict(
        #:frequency_threshold => 0,
        #:min_shape_r2 => 0.5,
        #:max_fwhm_ppm => 30.0,
        #:min_fwhm_ppm => 2.0
    #),
    #:PeakBinning => Dict(
        #:tolerance => 30.0,
        #:frequency_threshold => 0
    #)
)

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

function ensure_output_dir()
    if !isdir(OUTPUT_DIR)
        mkdir(OUTPUT_DIR)
        println("Created output directory: $OUTPUT_DIR")
    end
end

function plot_spectrum_step(mz, intensity, step_name; peaks=nothing, spectrum_index=1)
    fig = Figure(size=(1200, 600))
    ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Intensity", title="Step: $step_name (Spectrum $spectrum_index)")
    
    lines!(ax, mz, intensity, color=:blue, linewidth=1)
    
    if peaks !== nothing && !isempty(peaks)
        peak_mz = [p.mz for p in peaks]
        peak_intensity = [p.intensity for p in peaks]
        scatter!(ax, peak_mz, peak_intensity, color=:red, markersize=5)
    end
    
    filename = joinpath(OUTPUT_DIR, "spectrum_$(step_name)_index_$spectrum_index.png")
    save(filename, fig)
    println("  - Saved plot: $(basename(filename))")
end

function print_step_header(step_name::String)
    println("\n" * ">"^50)
    println("APPLYING STEP: $step_name")
    println(">"^50)
end

# ===================================================================
# PREPROCESSING PIPELINE FUNCTIONS (IN-PLACE)
# ===================================================================

function apply_baseline_correction_core(spectra::Vector{MutableSpectrum}, params::Dict, msi_data::MSIData)
    print_step_header("Baseline Correction")
    
    method = get(params, :method, :snip)
    iterations = get(params, :iterations, 100)
    window = get(params, :window, 20)
    println("  Method: $method, Iterations: $iterations, Window: $window")

    # Safely plot the first spectrum
    if !isempty(spectra)
        s = spectra[1]
        if validate_spectrum(s.mz, s.intensity)
            baseline = MSI_src.apply_baseline_correction_core(s.intensity; method=method, iterations=iterations, window=window)
            corrected_intensity = max.(0.0, s.intensity .- baseline)
            plot_spectrum_step(s.mz, corrected_intensity, "baseline_correction")
        end
    end

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            baseline = MSI_src.apply_baseline_correction_core(s.intensity; method=method, iterations=iterations, window=window)
            s.intensity = max.(0.0, s.intensity .- baseline)
        end
    end
end

function apply_smoothing(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Smoothing")
    
    method = get(params, :method, :savitzky_golay)
    window = get(params, :window, 9)
    order = get(params, :order, 2)
    println("  Method: $method, Window: $window, Order: $order")

    # Safely plot the first spectrum
    if !isempty(spectra)
        s = spectra[1]
        if validate_spectrum(s.mz, s.intensity)
            smoothed_intensity = max.(0.0, smooth_spectrum_core(s.intensity; method=method, window=window, order=order))
            plot_spectrum_step(s.mz, smoothed_intensity, "smoothing", spectrum_index=s.id)
        end
    end

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            smoothed_intensity = max.(0.0, smooth_spectrum_core(s.intensity; method=method, window=window, order=order))
            s.intensity = smoothed_intensity
        end
    end
end

function apply_peak_picking(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Peak Picking")
    
    method = get(params, :method, :profile)
    snr_threshold = get(params, :snr_threshold, 3.0)
    half_window = get(params, :half_window, 10)
    min_peak_prominence = get(params, :min_peak_prominence, 0.1)
    merge_peaks_tolerance = get(params, :merge_peaks_tolerance, 0.002)
    println("  Method: $method, SNR: $snr_threshold, Half Window: $half_window")

    found_invalid_spectrum = Threads.Atomic{Bool}(false)

    Threads.@threads for s in spectra
        is_valid = validate_spectrum(s.mz, s.intensity)

        if !is_valid && !found_invalid_spectrum[]
            if Threads.atomic_xchg!(found_invalid_spectrum, true) == false
                # Debug printing for the first invalid spectrum found
            end
        end

        if is_valid
            if method == :profile
                s.peaks = detect_peaks_profile_core(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window, min_peak_prominence=min_peak_prominence, merge_peaks_tolerance=merge_peaks_tolerance)
            elseif method == :wavelet
                s.peaks = detect_peaks_wavelet(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window)
            elseif method == :centroid
                s.peaks = detect_peaks_centroid_core(s.mz, s.intensity; snr_threshold=snr_threshold)
            else
                s.peaks = detect_peaks_profile_core(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window)
            end
        else
            s.peaks = []
        end
    end

    # Safely plot the first spectrum
    if !isempty(spectra)
        s1 = findfirst(s -> s.id == 1, spectra)
        if s1 !== nothing
            plot_spectrum_step(spectra[s1].mz, spectra[s1].intensity, "peak_picking", peaks=spectra[s1].peaks)
            println("  - Detected $(length(spectra[s1].peaks)) peaks in spectrum 1")
        end
    end
end

function apply_peak_selection(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Peak Selection")
    
    min_snr = get(params, :min_snr, 0.0)
    min_fwhm = get(params, :min_fwhm_ppm, 0.0)
    max_fwhm = get(params, :max_fwhm_ppm, Inf)
    min_r2 = get(params, :min_shape_r2, 0.0)
    
    # Handle `nothing` values from params, default to non-filtering values
    min_snr = isnothing(min_snr) ? 0.0 : min_snr
    min_fwhm = isnothing(min_fwhm) ? 0.0 : min_fwhm
    max_fwhm = isnothing(max_fwhm) ? Inf : max_fwhm
    min_r2 = isnothing(min_r2) ? 0.0 : min_r2

    println("  - Min SNR: $min_snr")
    println("  - FWHM Range (ppm): [$min_fwhm, $max_fwhm]")
    println("  - Min Shape R²: $min_r2")

    total_peaks_before = sum(s -> length(s.peaks), spectra)

    Threads.@threads for s in spectra
        if !isempty(s.peaks)
            s.peaks = filter(p -> 
                p.snr >= min_snr &&
                (min_fwhm <= p.fwhm <= max_fwhm) &&
                p.shape_r2 >= min_r2,
                s.peaks
            )
        end
    end

    total_peaks_after = sum(s -> length(s.peaks), spectra)
    println("  - Peaks before: $total_peaks_before, Peaks after: $total_peaks_after")

    # Safely plot the first spectrum to show effect of filtering
    if !isempty(spectra)
        s1_idx = findfirst(s -> s.id == 1, spectra)
        if s1_idx !== nothing
            s1 = spectra[s1_idx]
            plot_spectrum_step(s1.mz, s1.intensity, "peak_selection", peaks=s1.peaks, spectrum_index=s1.id)
            println("  - Spectrum 1 now has $(length(s1.peaks)) peaks after selection.")
        end
    end
end

function apply_calibration(spectra::Vector{MutableSpectrum}, params::Dict, reference_peaks::Dict)
    print_step_header("Calibration")
    
    method = get(params, :method, :none)
    ppm_tolerance = get(params, :ppm_tolerance, 20.0)
    println("  Method: $method, PPM Tolerance: $ppm_tolerance")

    if method == :none || isempty(reference_peaks)
        println("  - Skipping calibration (no method or reference peaks)")
        return
    end

    reference_masses = collect(keys(reference_peaks))
    calibration_info = Vector{String}(undef, length(spectra))

    Threads.@threads for i in 1:length(spectra)
        s = spectra[i]
        info_message = ""
        if validate_spectrum(s.mz, s.intensity)
            matched_peaks = find_calibration_peaks_core(s.mz, s.intensity, reference_masses; ppm_tolerance=ppm_tolerance)
            if length(matched_peaks) >= 2
                measured = sort(collect(values(matched_peaks)))
                theoretical = sort(collect(keys(matched_peaks)))
                itp = linear_interpolation(measured, theoretical, extrapolation_bc=Line())
                s.mz = itp(s.mz) # Modify mz-axis in-place
                info_message = "  - Spectrum $(s.id): calibrated using $(length(matched_peaks)) reference peaks"
            else
                info_message = "  - Spectrum $(s.id): insufficient reference peaks ($(length(matched_peaks)) found), skipping"
            end
        end
        calibration_info[i] = info_message
    end

    # Print summary of results
    printed_count = 0
    for info in calibration_info
        if !isempty(info) && printed_count < 5
            println(info)
            printed_count += 1
        end
    end
end

function apply_peak_alignment(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Peak Alignment")
    
    method = get(params, :method, :none)
    tolerance = get(params, :tolerance, 0.002)
    tolerance_unit = get(params, :tolerance_unit, :mz)
    println("  Method: $method, Tolerance: $tolerance $tolerance_unit")

    if method == :none
        println("  - Skipping peak alignment")
        return
    end

    ref_find_idx = findfirst(s -> !isempty(s.peaks), spectra)
    if ref_find_idx === nothing
        println("  - Insufficient spectra with peaks for alignment. Skipping.")
        return
    end

    ref_spectrum = spectra[ref_find_idx]
    ref_peaks_mz = [p.mz for p in ref_spectrum.peaks]
    println("  - Using spectrum $(ref_spectrum.id) as reference with $(length(ref_peaks_mz)) peaks")

    Threads.@threads for s in spectra
        if s.id == ref_spectrum.id || isempty(s.peaks)
            continue
        end
        
        current_peaks_mz = [p.mz for p in s.peaks]
        alignment_func = align_peaks_lowess_core(ref_peaks_mz, current_peaks_mz; method=method, tolerance=tolerance, tolerance_unit=tolerance_unit)
        
        s.mz = alignment_func.(s.mz) # Update m/z axis
        
        # Update peak m/z values
        for i in 1:length(s.peaks)
            old_peak = s.peaks[i]
            aligned_peak_mz = alignment_func(old_peak.mz)
            s.peaks[i] = (mz=aligned_peak_mz, intensity=old_peak.intensity, fwhm=old_peak.fwhm, shape_r2=old_peak.shape_r2, snr=old_peak.snr, prominence=old_peak.prominence)
        end
    end
end

function apply_normalization_core(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Normalization")
    
    method = get(params, :method, :tic)
    println("  Method: $method")

    # Safely plot the first spectrum
    if !isempty(spectra)
        s1 = findfirst(s -> s.id == 1, spectra)
        if s1 !== nothing
            s = spectra[s1]
            if validate_spectrum(s.mz, s.intensity)
                normalized_intensity = MSI_src.apply_normalization_core(s.intensity; method=method)
                plot_spectrum_step(s.mz, normalized_intensity, "normalization")
            end
        end
    end

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            s.intensity = MSI_src.apply_normalization_core(s.intensity; method=method)
        end
    end
end

function apply_peak_binning(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Peak Binning")
    
    method = get(params, :method, :adaptive)
    tolerance = get(params, :tolerance, 20.0)
    tolerance_unit = get(params, :tolerance_unit, :ppm)
    min_peak_per_bin = get(params, :min_peak_per_bin, 3)
    println("  Method: $method, Tolerance: $tolerance $tolerance_unit")

    if isempty(spectra) || all(s -> isempty(s.peaks), spectra)
        @warn "No peaks found for binning. Returning empty feature matrix."
        return nothing, nothing
    end
    
    println("  - Binning peaks from $(length(spectra)) spectra")
    
    # Collect all peaks with their intensities
    all_peaks = Vector{Tuple{Float64, Float64}}()  # (mz, intensity)
    for s in spectra
        for p in s.peaks
            push!(all_peaks, (p.mz, p.intensity))
        end
    end
    
    if isempty(all_peaks)
        @warn "No peaks collected for binning."
        return nothing, nothing
    end
    
    sort!(all_peaks, by=x->x[1])
    
    # Create bins - just store mz_center and intensity
    bin_centers = Float64[]
    bin_intensities = Float64[]
    
    i = 1
    while i <= length(all_peaks)
        current_bin_start = i
        current_peak = all_peaks[i]
        
        # Find all peaks in this bin
        j = i + 1
        while j <= length(all_peaks)
            next_peak = all_peaks[j]
            # Calculate tolerance
            tol = (tolerance_unit == :ppm) ? (current_peak[1] * tolerance / 1e6) : tolerance
            
            if (next_peak[1] - current_peak[1]) <= tol
                j += 1
            else
                break
            end
        end
        
        current_bin_end = j - 1
        bin_size = current_bin_end - current_bin_start + 1
        
        # Check if we have enough peaks in this bin
        if bin_size >= min_peak_per_bin
            bin_peaks_core = all_peaks[current_bin_start:current_bin_end]
            
            # Calculate m/z center and average intensity
            mz_sum = 0.0
            intensity_sum = 0.0
            for peak in bin_peaks_core
                mz_sum += peak[1]
                intensity_sum += peak[2]
            end
            
            mz_center = mz_sum / bin_size
            avg_intensity = intensity_sum / bin_size
            
            push!(bin_centers, mz_center)
            push!(bin_intensities, avg_intensity)
        end
        
        i = j  # Move to next potential bin
    end
    
    # Create the 2-row matrix
    if !isempty(bin_centers)
        n_bins = length(bin_centers)
        feature_matrix = Matrix{Float64}(undef, 2, n_bins)
        
        for i in 1:n_bins
            feature_matrix[1, i] = bin_centers[i]
            feature_matrix[2, i] = bin_intensities[i]
        end
        
        println("  - Created feature matrix: 2 × $n_bins")
        println("  - Number of bins created: $n_bins")
        println("  - m/z range: $(round(bin_centers[1], digits=4)) - $(round(bin_centers[end], digits=4))")
        
        # Return bin info as a vector of tuples for compatibility
        bin_info = [(bin_centers[i], bin_intensities[i]) for i in 1:n_bins]
        return feature_matrix, bin_info
    else
        @warn "No bins created after filtering"
        return nothing, nothing
    end
end

function save_feature_matrix(feature_matrix::Matrix{Float64}, bin_info)
    print_step_header("Saving Results")
    
    # Save as simple CSV with m/z and intensity rows
    csv_path = joinpath(OUTPUT_DIR, "feature_matrix_simple.csv")
    
    open(csv_path, "w") do io
        # Write header
        write(io, "mz,intensity\n")
        
        # Write data: m/z values in first column, intensities in second
        for i in 1:size(feature_matrix, 2)
            mz = feature_matrix[1, i]
            intensity = feature_matrix[2, i]
            write(io, "$mz,$intensity\n")
        end
    end
    println("  - Saved simple feature matrix: $csv_path")
    
    # Also save in a more standard format for MSI
    csv_path_standard = joinpath(OUTPUT_DIR, "feature_matrix_standard.csv")
    
    open(csv_path_standard, "w") do io
        # Write header with m/z values as column names
        write(io, "sample_type,")
        mz_headers = [@sprintf("mz_%.4f", feature_matrix[1, i]) for i in 1:size(feature_matrix, 2)]
        write(io, join(mz_headers, ",") * "\n")
        
        # Write the aggregated intensity values
        write(io, "aggregated_spectrum,")
        intensity_values = [feature_matrix[2, i] for i in 1:size(feature_matrix, 2)]
        write(io, join(string.(intensity_values), ",") * "\n")
    end
    println("  - Saved standard format matrix: $csv_path_standard")
    
    return csv_path, csv_path_standard
end

# ===================================================================
# MAIN PREPROCESSING PIPELINE
# ===================================================================

function run_preprocessing_pipeline()
    println("="^80)
    println("MSI PREPROCESSING PIPELINE")
    println("="^80)
    
    # Ensure output directory exists
    ensure_output_dir()
    
    # Load data
    println("\nLoading data: $(basename(TEST_FILE))")
    msi_data = OpenMSIData(TEST_FILE)

    # Precompute analytics
    println("\nPrecomputing analytics...")
    precompute_analytics(msi_data)
    
    # Run precalculation to get automatic parameters
    println("\nRunning precalculation for automatic parameter determination...")
    auto_params = main_precalculation(msi_data, reference_peaks=reference_peaks, mask_path=MASK_ROUTE)
    
    # --- Apply User Overrides ---
    if !isempty(USER_OVERRIDES)
        println("\nApplying user overrides...")
        for (step, params) in USER_OVERRIDES
            if haskey(auto_params, step)
                for (param, value) in params
                    @info "OVERRIDE: For step `:$step`, setting `:$param` to `$value` (was `$(get(auto_params[step], param, "not set"))`)"
                end
                merge!(auto_params[step], params)
            else
                @warn "User override for non-existent step `:$step` ignored."
            end
        end
    end

    # --- Final Parameter Summary ---
    println("\n" * "-"^60)
    println("FINAL PARAMETERS FOR PIPELINE RUN")
    println("-"^60)
    for (step, params) in sort(collect(pairs(auto_params)), by=p -> p.first)
        println("  - Step: $step")
        if isempty(params)
            println("    (No parameters)")
            continue
        end
        for (param, value) in sort(collect(pairs(params)), by=p -> p.first)
            println("    - $(rpad(param, 25)): $value")
        end
    end
    
    # Define pipeline steps
    pipeline_steps = PIPELINE_STP
    
    println("\nPipeline steps: $(join(pipeline_steps, " -> "))")
    
    # Initialize `current_spectra` as a Vector of mutable structs for in-place modification
    println("\nInitializing spectra data structure...")
    local spectrum_indices_to_process::AbstractVector{Int}
    if !isempty(MASK_ROUTE)
        println("Applying mask from: $(MASK_ROUTE)")
        try
            mask_matrix = MSI_src.load_and_prepare_mask(MASK_ROUTE, msi_data.image_dims)
            masked_indices_set = MSI_src.get_masked_spectrum_indices(msi_data, mask_matrix)
            spectrum_indices_to_process = collect(masked_indices_set)
            println("Mask applied. $(length(spectrum_indices_to_process)) spectra are within the masked region.")
        catch e
            @error "Failed to load or apply mask: $e. Proceeding without mask."
            spectrum_indices_to_process = 1:length(msi_data.spectra_metadata)
        end
    else
        spectrum_indices_to_process = 1:length(msi_data.spectra_metadata)
    end

    if isempty(spectrum_indices_to_process)
        @warn "No spectra available for processing after applying mask/filter. Exiting pipeline."
        close(msi_data)
        return
    end

    num_spectra_to_process = length(spectrum_indices_to_process)
    current_spectra = Vector{MutableSpectrum}(undef, num_spectra_to_process)
    
    Threads.@threads for i in 1:num_spectra_to_process
        original_idx = spectrum_indices_to_process[i]
        mz, intensity = MSI_src.GetSpectrum(msi_data, original_idx)
        current_spectra[i] = MSI_src.MutableSpectrum(original_idx, mz, intensity, [])
        
        # Plot raw spectrum without masks (only the first processed one)
        if i == 1
            plot_spectrum_step(mz, intensity, "raw_unmasked_spectrum", spectrum_index=original_idx)
        end
    end
    
    # These variables will be populated by the pipeline steps
    feature_matrix = nothing
    bin_definitions = nothing

    # Apply pipeline steps, modifying `current_spectra` in-place
    for step in pipeline_steps
        println("\n" * "-"^60)
        println("PROCESSING STEP: $step")
        println("-"^60)
        
        if step == "stabilization"
            print_step_header("Intensity Transformation (Stabilization)")
            method = get(auto_params[:Stabilization], :method, :sqrt)
            println("  Method: $method")
            
            # Safely plot the first spectrum before transformation
            if !isempty(current_spectra)
                s = current_spectra[1]
                if validate_spectrum(s.mz, s.intensity)
                    # Use a temporary spectrum to plot before and after
                    initial_intensity = deepcopy(s.intensity)
                    plot_spectrum_step(s.mz, initial_intensity, "stabilization_before", spectrum_index=s.id)
                end
            end

            @time apply_intensity_transformation(current_spectra, auto_params[:Stabilization])
            
            # Safely plot the first spectrum after transformation
            if !isempty(current_spectra)
                s = current_spectra[1]
                if validate_spectrum(s.mz, s.intensity)
                    plot_spectrum_step(s.mz, s.intensity, "stabilization_after", spectrum_index=s.id)
                end
            end
            
        elseif step == "baseline_correction"
            @time apply_baseline_correction_core(current_spectra, auto_params[:BaselineCorrection], msi_data)
            
        elseif step == "smoothing"
            apply_smoothing(current_spectra, auto_params[:Smoothing])
            
        elseif step == "peak_picking"
            @time apply_peak_picking(current_spectra, auto_params[:PeakPicking])

        elseif step == "peak_selection"
            @time apply_peak_selection(current_spectra, auto_params[:PeakSelection])
            
        elseif step == "calibration"
            @time apply_calibration(current_spectra, auto_params[:Calibration], reference_peaks)
            
        elseif step == "peak_alignment"
            @time apply_peak_alignment(current_spectra, auto_params[:PeakAlignment])
            
        elseif step == "normalization"
            @time apply_normalization_core(current_spectra, auto_params[:Normalization])
            
        elseif step == "peak_binning"
            # This step is different as it generates the final matrix, not modifying spectra in-place
            feature_matrix, bin_info = @time apply_peak_binning(current_spectra, auto_params[:PeakBinning])
            
            if feature_matrix !== nothing
                @time save_feature_matrix(feature_matrix, bin_info)
            end
            
        else
            @warn "Unknown step: $step, skipping"
        end
        
        # Print progress
        println("✓ Completed step: $step")
    end
    
    # Clean up
    close(msi_data)
    GC.gc()
    
    println("\n" * "="^80)
    println("PREPROCESSING PIPELINE COMPLETED SUCCESSFULLY!")
    println("Results saved to: $OUTPUT_DIR")
    println("="^80)
end

# ===================================================================
# EXECUTE PIPELINE
# ===================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    @time run_preprocessing_pipeline()
end