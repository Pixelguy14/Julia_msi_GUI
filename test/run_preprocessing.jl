# test/run_preprocessing.jl

# ===================================================================
# Preprocessing Test Suite for JuliaMSI
# ===================================================================
# This script provides a customizable workflow to test and visualize
# the effects of different preprocessing steps on individual spectra
# from .mzML or .imzML files.
#
# Instructions:
# 1. Configure the file paths and parameters in the "CONFIG" section.
# 2. Run the script from the project's root directory:
#    julia --threads auto --project=. test/run_preprocessing.jl
# 3. Check the configured `RESULTS_DIR` for output plots and CSV files.
# ===================================================================

using Printf
using CairoMakie
using DataFrames
using CSV
using Statistics
using Interpolations
import Pkg

# --- Load the MSI_src Module ---
# Activate the project environment at the parent directory of this test script
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src


# ===================================================================
# CONFIG: CUSTOMIZE YOUR PREPROCESSING WORKFLOW HERE
# ===================================================================

# --- Input and Output ---
const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML" # Can be .mzML or .imzML
const RESULTS_DIR = "test/results/preprocessing"
const NUM_SPECTRA_TO_PROCESS = nothing # Set to `nothing` to process all spectra

# --- Internal Standards for Calibration and QC ---
# replace these with m/z values of known compounds present in your dataset.
const INTERNAL_STANDARDS = Dict(
    "P13" => 432.6584,
    "P15" => 464.6059,
    "P17" => 526.5534,
    "P21" => 650.4485,
    "P25" => 774.3435,
    "P29" => 898.2385,
    "P31" => 950.1861,
    "P33" => 1022.1336,
    "P37" => 1146.0286,
    "P45" => 1593.8187,
    "Unknown 1" => 772.433,
    "Unknown 2" => 772.5253
)

# --- Preprocessing Step Parameters ---

# Step 0: Quality Control (QC)
const QC_PARAMS = (
    ppm_tolerance = 5.0, # PPM tolerance for matching internal standards
)

# Step 1: Calibration
const CALIBRATION_PARAMS = (
    enabled = true,
    ppm_tolerance = 5.0, # PPM tolerance for finding calibration peaks
)

# Step 2: Smoothing
# Note: Smoothing is less effective and often unnecessary for centroid data.
const SMOOTHING_PARAMS = (
    enabled = false,
    window = 9,
    order = 2,
)

# Step 3: Baseline Correction
# Note: Baseline correction is less effective and often unnecessary for centroid data.
const BASELINE_CORRECTION_PARAMS = (
    enabled = false,
    iterations = 100,
)

# Step 4: Normalization
const NORMALIZATION_PARAMS = (
    enabled = true,
    method = :tic, # :tic, :median, or :none
)

# Step 5: Peak Detection
const PEAK_DETECTION_PARAMS = (
    enabled = true,
    method = :centroid, # Use :centroid for centroided data, :profile for profile data
    # --- Parameters for :centroid method ---
    intensity_threshold = 0.0, # Filters out peaks below this absolute intensity
    # --- Parameters for :profile method ---
    half_window = 10,
    snr_threshold = 3.0,
    min_peak_prominence = 0.05,
    merge_peaks_tolerance = 0.002,
)

# Step 6: Peak Alignment (Warping)
# This is performed after collecting peaks from all spectra.
const PEAK_ALIGNMENT_PARAMS = (
    enabled = true,
    tolerance = 0.002,
    tolerance_unit = :mz, # :mz or :ppm
    min_matched_peaks = 3, # Lowered for sparse centroid data
)

# Step 7: Peak Binning
const PEAK_BINNING_PARAMS = (
    enabled = true,
    tolerance = 0.1,
    tolerance_unit = :mz, # :mz or :ppm
    frequency_threshold = 0.1, # Min fraction of spectra a peak must be in
)

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

"""
    plot_spectrum_on_fig(fig_pos, mz, intensity, title)

Helper function to plot a spectrum on a specific position of a CairoMakie figure.
"""
function plot_spectrum_on_fig(fig_pos, mz, intensity, title)
    ax = Axis(fig_pos, title=title, xlabel="m/z", ylabel="Intensity")
    lines!(ax, mz, intensity)
end


# ===================================================================
# MAIN PROCESSING SCRIPT
# ===================================================================

function run_preprocessing_suite()
    println("Starting Preprocessing Test Suite...")
    mkpath(RESULTS_DIR)

    # --- Load Data ---
    println("Loading data from: $TEST_FILE")
    if !isfile(TEST_FILE)
        println("ERROR: Test file not found. Please check the path in the CONFIG section.")
        return
    end
    msi_data = @time OpenMSIData(TEST_FILE)
    println("Data loaded successfully. Found $(length(msi_data.spectra_metadata)) spectra.")

    # --- Determine which spectra to process ---
    total_spectra = length(msi_data.spectra_metadata)
    indices_to_process = if NUM_SPECTRA_TO_PROCESS === nothing
        1:total_spectra
    else
        unique(round.(Int, range(1, total_spectra, length=min(NUM_SPECTRA_TO_PROCESS, total_spectra))))
    end
    println("Will process $(length(indices_to_process)) spectra.")

    # --- Data storage for results ---
    qc_results = DataFrame(spectrum_idx=Int[], metric=String[], value=Float64[], compound=String[])
    all_processed_peaks_mz = Vector{Vector{Float64}}()
    all_processed_peaks_int = Vector{Vector{Float64}}()
    spectrum_indices_with_peaks = Int[]

    # --- Main Loop: Process each spectrum individually ---
    println("\n" * "="^20 * " Processing Individual Spectra " * "="^20)
    for (i, idx) in enumerate(indices_to_process)
        print("\rProcessing spectrum #$idx ($(i)/$(length(indices_to_process)))...")
        
        process_spectrum(msi_data, idx) do mz, intensity
            if qc_is_empty(mz, intensity) || !qc_is_regular(mz)
                # @warn "Skipping empty or irregular spectrum #$idx"
                return
            end

            original_mz, original_intensity = copy(mz), copy(intensity)
            processed_mz, processed_intensity = copy(mz), copy(intensity)

            # --- Step 0: QC (PPM and Resolution) ---
            # This happens before any modification to the m/z axis
            let ref_masses = collect(values(INTERNAL_STANDARDS))
                matched_peaks = find_calibration_peaks(processed_mz, processed_intensity, ref_masses, ppm_tolerance=QC_PARAMS.ppm_tolerance)
                
                for (compound, theoretical_mz) in INTERNAL_STANDARDS
                    # Find the measured m/z that corresponds to this theoretical_mz
                    measured_mz = 0.0
                    for (theo, meas) in matched_peaks
                        if theo == theoretical_mz
                            measured_mz = meas
                            break
                        end
                    end

                    if measured_mz > 0
                        # PPM Error
                        ppm_error = calculate_ppm_error(measured_mz, theoretical_mz)
                        push!(qc_results, (idx, "ppm_error", ppm_error, compound))
                        
                        # Resolution
                        resolution = calculate_resolution_fwhm(measured_mz, processed_mz, processed_intensity)
                        if !isnan(resolution)
                            push!(qc_results, (idx, "resolution", resolution, compound))
                        end
                    end
                end
            end

            # --- Step 1: Calibration ---
            if CALIBRATION_PARAMS.enabled
                ref_masses = collect(values(INTERNAL_STANDARDS))
                matched_peaks = find_calibration_peaks(processed_mz, processed_intensity, ref_masses, ppm_tolerance=CALIBRATION_PARAMS.ppm_tolerance)
                if length(matched_peaks) >= 2
                    measured = sort(collect(values(matched_peaks)))
                    theoretical = sort(collect(keys(matched_peaks)))
                    itp = linear_interpolation(measured, theoretical, extrapolation_bc=Line())
                    processed_mz = itp(processed_mz)
                end
            end

            # --- Step 2: Smoothing ---
            if SMOOTHING_PARAMS.enabled
                processed_intensity = smooth_spectrum(processed_intensity, window=SMOOTHING_PARAMS.window, order=SMOOTHING_PARAMS.order)
            end

            # --- Step 3: Baseline Correction ---
            if BASELINE_CORRECTION_PARAMS.enabled
                baseline = snip_baseline(processed_intensity, iterations=BASELINE_CORRECTION_PARAMS.iterations)
                processed_intensity .-= baseline
                processed_intensity = max.(0, processed_intensity) # Ensure non-negativity
            end

            # --- Step 4: Normalization ---
            if NORMALIZATION_PARAMS.enabled && NORMALIZATION_PARAMS.method != :none
                if NORMALIZATION_PARAMS.method == :tic
                    processed_intensity = tic_normalize(processed_intensity)
                elseif NORMALIZATION_PARAMS.method == :median
                    processed_intensity = median_normalize(processed_intensity)
                end
            end

            # --- Step 5: Peak Detection ---
            local pk_mz, pk_int
            pk_mz, pk_int = Float64[], Float64[] # Initialize empty
            if PEAK_DETECTION_PARAMS.enabled
                if PEAK_DETECTION_PARAMS.method == :profile
                    pk_mz, pk_int = detect_peaks_profile(processed_mz, processed_intensity,
                        half_window=PEAK_DETECTION_PARAMS.half_window,
                        snr_threshold=PEAK_DETECTION_PARAMS.snr_threshold,
                        min_peak_prominence=PEAK_DETECTION_PARAMS.min_peak_prominence,
                        merge_peaks_tolerance=PEAK_DETECTION_PARAMS.merge_peaks_tolerance)
                elseif PEAK_DETECTION_PARAMS.method == :wavelet
                    pk_mz, pk_int = detect_peaks_wavelet(processed_mz, processed_intensity)
                else # :centroid
                    pk_mz, pk_int = detect_peaks_centroid(processed_mz, processed_intensity)
                end
                
                if !isempty(pk_mz)
                    push!(all_processed_peaks_mz, pk_mz)
                    push!(all_processed_peaks_int, pk_int)
                    push!(spectrum_indices_with_peaks, idx)
                end
            end

            # --- Visualization of a sample spectrum ---
            if i == 1 # Only plot the first processed spectrum
                println("\nGenerating example plots for spectrum #$idx...")
                fig = Figure(size=(1200, 800))
                
                plot_spectrum_on_fig(fig[1,1], original_mz, original_intensity, "1. Original Spectrum")
                plot_spectrum_on_fig(fig[2,1], processed_mz, processed_intensity, "2. After All Steps (Before Peak Picking)")
                
                # Plot detected peaks
                ax = Axis(fig[3,1], title="3. Detected Peaks")
                if !isempty(pk_mz)
                    stem!(ax, pk_mz, pk_int)
                end
                
                save(joinpath(RESULTS_DIR, "example_spectrum_processing.png"), fig)
                println("Saved example processing plots.")
            end
        end
    end
    println("\nIndividual spectrum processing complete.")

    # --- Save QC Results ---
    if !isempty(qc_results)
        println("\n" * "="^20 * " Generating QC Report " * "="^20)
        CSV.write(joinpath(RESULTS_DIR, "qc_results.csv"), qc_results)
        println("QC results saved to qc_results.csv")

        # Generate summary plots for QC
        fig = Figure(size=(1200, 600))
        
        # PPM Error Histogram
        ppm_errors = filter(row -> row.metric == "ppm_error", qc_results).value
        if !isempty(ppm_errors)
            ax1 = Axis(fig[1,1], title="PPM Error Distribution", xlabel="PPM Error")
            hist!(ax1, ppm_errors, bins=30)
        end

        # Resolution Histogram
        resolutions = filter(row -> row.metric == "resolution", qc_results).value
        if !isempty(resolutions)
            ax2 = Axis(fig[1,2], title="Resolution Distribution", xlabel="Resolution (FWHM)")
            hist!(ax2, resolutions, bins=30)
        end
        
        save(joinpath(RESULTS_DIR, "qc_summary_plots.png"), fig)
        println("QC summary plots saved.")
    end

    if isempty(all_processed_peaks_mz)
        @warn "No peaks were detected in any of the processed spectra. Skipping alignment and binning."
        return
    end

    # --- Step 6: Peak Alignment ---
    println("\n" * "="^20 * " Aligning Peaks " * "="^20)
    aligned_peaks_mz = copy(all_processed_peaks_mz)
    if PEAK_ALIGNMENT_PARAMS.enabled
        # Create a reference peak list (e.g., from the spectrum with the most peaks)
        ref_idx = argmax(length.(all_processed_peaks_mz))
        reference_peaks = all_processed_peaks_mz[ref_idx]
        println("Using spectrum $(spectrum_indices_with_peaks[ref_idx]) as alignment reference.")

        for i in 1:length(aligned_peaks_mz)
            if i == ref_idx continue end
            alignment_func = align_peaks_lowess(reference_peaks, aligned_peaks_mz[i],
                tolerance=PEAK_ALIGNMENT_PARAMS.tolerance,
                tolerance_unit=PEAK_ALIGNMENT_PARAMS.tolerance_unit,
                min_matched_peaks=PEAK_ALIGNMENT_PARAMS.min_matched_peaks)
            
            aligned_peaks_mz[i] = alignment_func(aligned_peaks_mz[i])
        end
        println("Peak alignment complete.")
    end

    # --- Step 7: Peak Binning & Feature Matrix Generation ---
    println("\n" * "="^20 * " Binning Peaks " * "="^20)
    if PEAK_BINNING_PARAMS.enabled
        feature_matrix, mz_bins = bin_peaks(aligned_peaks_mz, all_processed_peaks_int,
            PEAK_BINNING_PARAMS.tolerance,
            tolerance_unit=PEAK_BINNING_PARAMS.tolerance_unit,
            frequency_threshold=PEAK_BINNING_PARAMS.frequency_threshold)
        
        if !isempty(feature_matrix)
            df = DataFrame(feature_matrix, :auto)
            rename!(df, ["bin_$(i)" for i in 1:size(df, 2)])
            insertcols!(df, 1, :spectrum_idx => spectrum_indices_with_peaks)
            
            CSV.write(joinpath(RESULTS_DIR, "feature_matrix.csv"), df)
            println("Feature matrix saved to feature_matrix.csv")
            
            # Save bin m/z ranges
            bin_df = DataFrame(bin_index=1:length(mz_bins), mz_start=[b[1] for b in mz_bins], mz_end=[b[2] for b in mz_bins])
            CSV.write(joinpath(RESULTS_DIR, "bin_definitions.csv"), bin_df)
            println("Bin definitions saved to bin_definitions.csv")
        else
            @warn "Feature matrix was empty after binning."
        end
    end

    println("\nPreprocessing Test Suite finished successfully!")
end


# --- Execute ---
run_preprocessing_suite()
