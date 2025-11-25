# run_preprocessing.jl

using Printf
import Pkg
using CairoMakie
using DataFrames    # For creating dataframes
using CSV

# --- Load the MSI_src Module ---
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src


# ===================================================================
# CONFIG: PLEASE FILL IN YOUR FILE PATHS HERE
# ===================================================================

const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Thricoderma_etc/Imaging_interaccion_trichoderma_vs_streptomyces.imzML"
# const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/set de datos MS/Atropina_tuneo_fraq_20ev.mzML"
# const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML"

# const MASK_ROUTE = "/home/pixel/Documents/Cinvestav_2025/JuliaMSI/public/css/masks/Stomach_DHB_uncompressed.png"
const MASK_ROUTE = "" 

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

function apply_baseline_correction(spectra::Vector{MutableSpectrum}, params::Dict, msi_data::MSIData)
    print_step_header("Baseline Correction")
    
    method = get(params, :method, :snip)
    iterations = get(params, :iterations, 100)
    window = get(params, :window, 20)
    println("  Method: $method, Iterations: $iterations, Window: $window")

    # Safely plot the first spectrum
    if !isempty(spectra)
        s = spectra[1]
        if validate_spectrum(s.mz, s.intensity)
            baseline = MSI_src.apply_baseline_correction(s.intensity; method=method, iterations=iterations, window=window)
            corrected_intensity = max.(0.0, s.intensity .- baseline)
            plot_spectrum_step(s.mz, corrected_intensity, "baseline_correction")
        end
    end

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            baseline = MSI_src.apply_baseline_correction(s.intensity; method=method, iterations=iterations, window=window)
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
            smoothed_intensity = max.(0.0, smooth_spectrum(s.intensity; method=method, window=window, order=order))
            plot_spectrum_step(s.mz, smoothed_intensity, "smoothing", spectrum_index=s.id)
        end
    end

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            smoothed_intensity = max.(0.0, smooth_spectrum(s.intensity; method=method, window=window, order=order))
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
                s.peaks = detect_peaks_profile(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window, min_peak_prominence=min_peak_prominence, merge_peaks_tolerance=merge_peaks_tolerance)
            elseif method == :wavelet
                s.peaks = detect_peaks_wavelet(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window)
            elseif method == :centroid
                s.peaks = detect_peaks_centroid(s.mz, s.intensity; snr_threshold=snr_threshold)
            else
                s.peaks = detect_peaks_profile(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window)
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
            matched_peaks = find_calibration_peaks(s.mz, s.intensity, reference_masses; ppm_tolerance=ppm_tolerance)
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
        alignment_func = align_peaks_lowess(ref_peaks_mz, current_peaks_mz; method=method, tolerance=tolerance, tolerance_unit=tolerance_unit)
        
        s.mz = alignment_func.(s.mz) # Update m/z axis
        
        # Update peak m/z values
        for i in 1:length(s.peaks)
            old_peak = s.peaks[i]
            aligned_peak_mz = alignment_func(old_peak.mz)
            s.peaks[i] = (mz=aligned_peak_mz, intensity=old_peak.intensity, fwhm=old_peak.fwhm, shape_r2=old_peak.shape_r2, snr=old_peak.snr, prominence=old_peak.prominence)
        end
    end
end

function apply_normalization(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Normalization")
    
    method = get(params, :method, :tic)
    println("  Method: $method")

    # Safely plot the first spectrum
    if !isempty(spectra)
        s1 = findfirst(s -> s.id == 1, spectra)
        if s1 !== nothing
            s = spectra[s1]
            if validate_spectrum(s.mz, s.intensity)
                normalized_intensity = MSI_src.apply_normalization(s.intensity; method=method)
                plot_spectrum_step(s.mz, normalized_intensity, "normalization")
            end
        end
    end

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            s.intensity = MSI_src.apply_normalization(s.intensity; method=method)
        end
    end
end

function apply_peak_binning(spectra::Vector{MutableSpectrum}, params::Dict)
    print_step_header("Peak Binning")
    
    method = get(params, :method, :adaptive)
    tolerance = get(params, :tolerance, 20.0)
    tolerance_unit = get(params, :tolerance_unit, :ppm)
    min_peak_per_bin = get(params, :min_peak_per_bin, 3)
    max_bin_width_ppm = get(params, :max_bin_width_ppm, 150.0)
    intensity_weighted_centers = get(params, :intensity_weighted_centers, true)
    println("  Method: $method, Tolerance: $tolerance $tolerance_unit")

    if isempty(spectra) || all(s -> isempty(s.peaks), spectra)
        @warn "No peaks found for binning. Returning empty feature matrix."
        return nothing, nothing
    end
    
    println("  - Binning peaks from $(length(spectra)) spectra")
    
    binning_params = PeakBinningParams(
        method=method,
        tolerance=tolerance,
        tolerance_unit=tolerance_unit,
        min_peak_per_bin=min_peak_per_bin,
        max_bin_width_ppm=max_bin_width_ppm,
        intensity_weighted_centers=intensity_weighted_centers
    )
    
    feature_matrix, bin_definitions = bin_peaks(spectra, binning_params)
    
    if feature_matrix !== nothing && bin_definitions !== nothing
        println("  - Generated feature matrix: $(size(feature_matrix.matrix))")
        println("  - Number of bins: $(length(bin_definitions))")
    end
    
    return feature_matrix, bin_definitions
end

function save_feature_matrix(feature_matrix, bin_definitions)
    print_step_header("Saving Results")
    
    # Save feature matrix as CSV
    csv_path = joinpath(OUTPUT_DIR, "feature_matrix.csv")
    
    # Create column headers
    bin_headers = ["bin_$(i)_$(round(def[1], digits=4))-$(round(def[2], digits=4))" 
                         for (i, def) in enumerate(bin_definitions)]
    
    # Open file for writing
    open(csv_path, "w") do io
        # Write header row
        write(io, "spectrum_index," * join(bin_headers, ",") * "\n")
        
        # Write data rows
        for r_idx in 1:size(feature_matrix.matrix, 1)
            write(io, "$(feature_matrix.sample_ids[r_idx]),")
            for c_idx in 1:size(feature_matrix.matrix, 2)
                write(io, "$(feature_matrix.matrix[r_idx, c_idx])")
                if c_idx < size(feature_matrix.matrix, 2)
                    write(io, ",")
                end
            end
            write(io, "\n")
        end
    end
    println("  - Saved feature matrix: $csv_path")
    
    # Save bin definitions (still using DataFrame as it's small)
    bins_path = joinpath(OUTPUT_DIR, "bin_definitions.csv")
    bins_df = DataFrame(
        bin_index = 1:length(bin_definitions),
        mz_start = [def[1] for def in bin_definitions],
        mz_end = [def[2] for def in bin_definitions],
        mz_center = [(def[1] + def[2])/2 for def in bin_definitions]
    )
    CSV.write(bins_path, bins_df)
    println("  - Saved bin definitions: $bins_path")
    
    return csv_path, bins_path
end

# ===================================================================
# USER OVERRIDES: Manually specify parameters here
# ===================================================================
# This dictionary allows you to override any auto-detected parameters.
# The structure should match the output of `main_precalculation`.
# Example: Force a less aggressive peak prominence threshold.
const USER_OVERRIDES = Dict(
    :PeakPicking => Dict(
        #:snr_threshold => 1.5,
        :half_window => 5,
        :snr_threshold => 15.0,
        #:merge_peaks_tolerance => 2.5,
        #:half_window => 2
    ),
    # :Smoothing => Dict(
    #     :window => 11
    # )
    #:PeakAlignment => Dict(
    #    :tolerance => 0.01
    #)
    :PeakBinningParams => Dict(
        :tolerance => 30.0
    )
)

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
    pipeline_steps = [
        "baseline_correction",
        "smoothing", 
        "peak_picking",
        "calibration",
        "peak_alignment",
        "normalization",
        "peak_binning"
    ]
    
    println("\nPipeline steps: $(join(pipeline_steps, " -> "))")
    
    # Initialize `current_spectra` as a Vector of mutable structs for in-place modification
    println("\nInitializing spectra data structure...")
    num_spectra = length(msi_data.spectra_metadata)
    current_spectra = Vector{MutableSpectrum}(undef, num_spectra)
    _iterate_spectra_fast(msi_data) do idx, mz, intensity
        current_spectra[idx] = MutableSpectrum(idx, mz, intensity, [])
        
        # Plot raw spectrum without masks
        if idx == 1
            plot_spectrum_step(mz, intensity, "raw_unmasked_spectrum")
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
        
        if step == "baseline_correction"
            @time apply_baseline_correction(current_spectra, auto_params[:BaselineCorrection], msi_data)
            
        elseif step == "smoothing"
            apply_smoothing(current_spectra, auto_params[:Smoothing])
            
        elseif step == "peak_picking"
            @time apply_peak_picking(current_spectra, auto_params[:PeakPicking])
            
        elseif step == "calibration"
            @time apply_calibration(current_spectra, auto_params[:Calibration], reference_peaks)
            
        elseif step == "peak_alignment"
            @time apply_peak_alignment(current_spectra, auto_params[:PeakAlignment])
            
        elseif step == "normalization"
            @time apply_normalization(current_spectra, auto_params[:Normalization])
            
        elseif step == "peak_binning"
            # This step is different as it generates the final matrix, not modifying spectra in-place
            feature_matrix, bin_definitions = @time apply_peak_binning(current_spectra, auto_params[:PeakBinningParams])
            
            if feature_matrix !== nothing
                @time save_feature_matrix(feature_matrix, bin_definitions)
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