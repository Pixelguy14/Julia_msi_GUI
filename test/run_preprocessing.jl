# test/run_preprocessing.jl

# ===================================================================
# Test Environment for the Preprocessing.jl Module
# ===================================================================
# This script tests the full preprocessing pipeline on single spectra
# and total spectra from both .mzML and .imzML files.
# It generates an overlay plot showing all preprocessing stages and
# saves the resulting feature matrix to a CSV file.
#
# Instructions:
# 1. Ensure the file paths in the "CONFIG" section are correct.
# 2. Run the script from the project's root directory:
#    julia test/run_preprocessing.jl
# 3. Check the `test/results/` folder for output plots and CSVs.
# ===================================================================

using Printf
using CairoMakie
import Pkg
using DataFrames # For saving FeatureMatrix to CSV
using CSV        # For saving FeatureMatrix to CSV
using Statistics # For mean()

# --- Load Modules ---
# Activate the project environment to access dependencies
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src # This brings in Preprocessing.jl functions via export

# ===================================================================
# CONFIG: Test files and parameters
# ===================================================================

# --- Test Files ---
# An mzML file for testing spectrum-based processing
# const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.mzML"
# const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/set de datos MS/Leaf_profile_LD_LTP_MS.mzML"
const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/set de datos MS/Escopolamina_tuneo_fraq_20ev.mzML"
#const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/set de datos MS/Atropina_tuneo_fraq_20ev.mzML"

const MZML_SPECTRUM_ID = 1

# An imzML file for testing
# const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.imzML"
# const IMZML_COORDS = (50, 50) 
const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML"
const IMZML_COORDS = (1997, 639) 

# --- Output Directory ---
const RESULTS_DIR = "test/results"

# ===================================================================
# HELPER FUNCTIONS FOR PLOTTING
# ===================================================================

"""
    plot_overlay_stages(collected_data, output_path, title)

Creates a single plot overlaying spectra from different preprocessing stages.
"""
function plot_overlay_stages(collected_data, output_path, title)
    fig = Figure(size = (1400, 800))
    ax = Axis(fig[1, 1], title=title, xlabel="m/z", ylabel="Intensity")

    colors = Makie.wong_colors() # A good set of distinct colors

    for (i, (stage, mz, intensity)) in enumerate(collected_data)
        color = colors[mod1(i, length(colors))] # Cycle through colors
        
        # Plot the spectrum as a line
        lines!(ax, mz, intensity, color=color, label=string(stage))

        # If it's the peaks stage, also mark the peak tops
        if stage == :peaks
            scatter!(ax, mz, intensity, color=color, marker=:circle, markersize=8, label="$(string(stage)) (tops)")
        end
    end
    axislegend(ax, position=:rt) # Right top position
    save(output_path, fig)
    println("SUCCESS: Overlay plot saved to $output_path")
end


# ===================================================================
# TEST DEFINITIONS
# ===================================================================

"""
    test_full_pipeline(msi_data, spectrum_id; output_dir, file_type_prefix)

Tests the full preprocessing pipeline on a single spectrum and saves a plot
for each intermediate step using the `on_stage` callback.
`spectrum_id` can be an `Int` (for mzML) or a `Tuple{Int, Int}` (for imzML).
"""
function test_full_pipeline(msi_data, spectrum_id; output_dir, file_type_prefix, mz_tolerance=0.002)
    println("\n--- Testing Full Preprocessing Pipeline on Spectrum: $spectrum_id (File Type: $file_type_prefix) ---")

    # 1. Determine the spectrum index
    local spec_idx
    if spectrum_id isa Int
        spec_idx = spectrum_id
    else # Tuple for imzML
        spec_idx = msi_data.coordinate_map[spectrum_id...]
    end

    if spec_idx == 0
        println("SKIPPED: No spectrum found at coordinates $spectrum_id.")
        return
    end

    # 2. Define the pipeline steps in the desired order
    pipeline_steps = [
        :qc,
        :transform,
        :smooth,
        :baseline,
        :normalize,
        :peaks,
        :align, # Align requires multiple spectra, but we'll run it on a single one for now (will warn)
        :bin
    ]

    # Define parameters for each step
    params = Dict(
        :transform_method => :sqrt,
        :sg_window => 15,
        :sg_order => 2,
        :snip_iterations => 100,
        :normalize_method => :tic,
        :peak_half_window => 10,
        :peak_snr => 3.0,
        :peak_intensity_threshold => 0.0, # For centroid peak detection
        :align_tolerance => mz_tolerance,
        :bin_tolerance => mz_tolerance,
        :bin_min_frequency => 0.0 # Keep all bins for a single spectrum
    )

    # 3. Define the on_stage callback to collect data for overlay plot and save separate plots
    collected_stage_data = []
    stage_counter = Ref(0) # Initialize counter for sequential naming
    normalized_spectrum = nothing # Variable to hold the normalized spectrum

    function stage_callback(stage; idx, mz, intensity)
        stage_counter[] += 1 # Increment counter
        println("  -> Generating plot for stage: $stage")
        
        local fig # Make fig available in the whole function scope

        if stage == :normalize
            normalized_spectrum = (mz, intensity)
            fig = plot_stage_spectrum(mz, intensity, title="Stage: $stage (Spectrum $spectrum_id)")
        elseif stage == :peaks && normalized_spectrum !== nothing
            # For the peaks stage, plot the normalized spectrum as a base layer
            fig = Figure(size = (1400, 500))
            ax = Axis(fig[1, 1], title="Stage: Peaks (Spectrum $spectrum_id)", xlabel="m/z", ylabel="Intensity")
            lines!(ax, normalized_spectrum[1], normalized_spectrum[2], color=:gray, label="Normalized Spectrum")
            scatter!(ax, mz, intensity, color=:red, marker=:circle, markersize=8, label="Detected Peaks")
            axislegend(ax)
        else
            # Default plotting for all other stages
            fig = plot_stage_spectrum(mz, intensity, title="Stage: $stage (Spectrum $spectrum_id)")
        end

        # Save the figure
        stage_output_path = joinpath(output_dir, "$(file_type_prefix)_$(spectrum_id)_$(stage_counter[])_$(stage).png")
        save(stage_output_path, fig)

        # Collect data for overlay plot
        push!(collected_stage_data, (stage, mz, intensity))
    end

    # 4. Run the pipeline on the single spectrum
    println("Running pipeline with steps: $pipeline_steps")
    processed_result = run_preprocessing_pipeline(
        msi_data, 
        [spec_idx], # The pipeline expects a vector of indices
        steps=pipeline_steps,
        params=params,
        on_stage=stage_callback
    )

    # 5. Generate and save the overlay plot
    overlay_output_path = joinpath(output_dir, "$(file_type_prefix)_$(spectrum_id)_all_stages_overlay.png")
    plot_overlay_stages(collected_stage_data, overlay_output_path, "Preprocessing Stages Overlay (Spectrum $spectrum_id)")

    # 6. Save feature matrix if generated
    if processed_result isa FeatureMatrix
        feature_matrix_output_path = joinpath(output_dir, "$(file_type_prefix)_$(spectrum_id)_feature_matrix.csv")
        # Convert mz_bins to a more readable format for CSV
        mz_labels = ["$(round(b[1], digits=4))_$(round(b[2], digits=4))" for b in processed_result.mz_bins]
        df = DataFrame(processed_result.matrix, Symbol.(mz_labels))
        CSV.write(feature_matrix_output_path, df)
        println("SUCCESS: Feature matrix saved to $feature_matrix_output_path")
    else
        @warn "Pipeline did not return a FeatureMatrix for Spectrum $spectrum_id."
        processed_result
    end

    println("--- Pipeline test finished for Spectrum: $spectrum_id (File Type: $file_type_prefix) ---")
    println("Check the '$(output_dir)' directory for output plots and CSVs.")
end

"""
    test_full_pipeline_on_total_spectrum(msi_data; output_dir, file_type_prefix)

Tests the full preprocessing pipeline on the *total spectrum* (sum of all spectra)
and saves a plot for each intermediate step.
"""
function test_full_pipeline_on_total_spectrum(msi_data; output_dir, file_type_prefix, mz_tolerance=0.002)
    println("\n--- Testing Full Preprocessing Pipeline on TOTAL Spectrum (File Type: $file_type_prefix) ---")

    # 1. Get the total spectrum
    total_mz, total_intensity = get_total_spectrum(msi_data)
    total_spectrum = (total_mz, total_intensity)

    if qc_is_empty(total_mz, total_intensity)
        println("SKIPPED: Total spectrum is empty.")
        return
    end

    # 2. Define the pipeline steps and parameters (same as for single spectrum)
    pipeline_steps = [
        :qc,
        :transform,
        :smooth,
        :baseline,
        :normalize,
        :peaks,
        :align, # Align requires multiple spectra, but we'll run it on a single one for now (will warn)
        :bin
    ]

    params = Dict(
        :transform_method => :sqrt,
        :sg_window => 15,
        :sg_order => 2,
        :snip_iterations => 100,
        :normalize_method => :tic,
        :peak_half_window => 10,
        :peak_snr => 3.0,
        :align_tolerance => mz_tolerance,
        :bin_tolerance => mz_tolerance,
        :bin_min_frequency => 0.0 # Keep all bins for a single spectrum
    )

    # 3. Define the on_stage callback
    collected_stage_data = []
    stage_counter = Ref(0) # Initialize counter for sequential naming
    normalized_spectrum_total = nothing # Variable to hold the normalized spectrum

    function stage_callback_total(stage; idx, mz, intensity)
        stage_counter[] += 1 # Increment counter
        println("  -> Generating plot for stage: $stage (Total Spectrum)")
        
        local fig

        if stage == :normalize
            normalized_spectrum_total = (mz, intensity)
            fig = plot_stage_spectrum(mz, intensity, title="Stage: $stage (Total Spectrum)")
        elseif stage == :peaks && normalized_spectrum_total !== nothing
            fig = Figure(size = (1400, 500))
            ax = Axis(fig[1, 1], title="Stage: Peaks (Total Spectrum)", xlabel="m/z", ylabel="Intensity")
            lines!(ax, normalized_spectrum_total[1], normalized_spectrum_total[2], color=:gray, label="Normalized Spectrum")
            scatter!(ax, mz, intensity, color=:red, marker=:circle, markersize=8, label="Detected Peaks")
            axislegend(ax)
        else
            fig = plot_stage_spectrum(mz, intensity, title="Stage: $stage (Total Spectrum)")
        end

        # Save separate plot with sequential name
        stage_output_path = joinpath(output_dir, "$(file_type_prefix)_total_$(stage_counter[])_$(stage).png")
        save(stage_output_path, fig)

        # Collect data for overlay plot
        push!(collected_stage_data, (stage, mz, intensity))
    end

    # 4. Run the pipeline on the single total spectrum
    println("Running pipeline with steps: $pipeline_steps")
    processed_result = run_preprocessing_pipeline(
        [total_spectrum], # Pass the total spectrum as a vector of one spectrum
        steps=pipeline_steps,
        params=params,
        on_stage=stage_callback_total
    )

    # 5. Generate and save the overlay plot
    overlay_output_path = joinpath(output_dir, "$(file_type_prefix)_total_all_stages_overlay.png")
    plot_overlay_stages(collected_stage_data, overlay_output_path, "Preprocessing Stages Overlay (Total Spectrum)")

    # 6. Save feature matrix if generated
    if processed_result isa FeatureMatrix
        feature_matrix_output_path = joinpath(output_dir, "$(file_type_prefix)_total_feature_matrix.csv")
        # Convert mz_bins to a more readable format for CSV
        mz_labels = ["$(round(b[1], digits=4))_$(round(b[2], digits=4))" for b in processed_result.mz_bins]
        df = DataFrame(processed_result.matrix, Symbol.(mz_labels))
        CSV.write(feature_matrix_output_path, df)
        println("SUCCESS: Feature matrix saved to $feature_matrix_output_path")
    else
        @warn "Pipeline did not return a FeatureMatrix for Total Spectrum."
        processed_result
    end

    println("--- Pipeline test finished for TOTAL Spectrum (File Type: $file_type_prefix) ---")
    println("Check the '$(output_dir)' directory for output plots and CSVs.")
end


# ===================================================================
# TEST RUNNER
# ===================================================================

function run_preprocessing_tests()
    println("="^80)
    println("STARTING PREPROCESSING TEST SUITE")
    println("="^80)

    # --- Test Case 1: Run full pipeline on a single mzML spectrum ---
    println("\n" * "="^20 * " Test Case 1: Full Pipeline on .mzML Spectrum " * "="^20)
    println("FILE: ", TEST_MZML_FILE)
    if isfile(TEST_MZML_FILE)
        try
            msi_data_mzml = OpenMSIData(TEST_MZML_FILE)
            
            # Dynamically determine tolerance
            println("\n--- Calculating optimal tolerance for .mzML data ---")
            report_mzml = analyze_mass_accuracy(msi_data_mzml, get_common_calibration_standards(:maldi_pos))
            mz_tolerance_mzml = 0.002 # Default
            if haskey(report_mzml, :optimal_ppm) && !isnan(report_mzml.optimal_ppm) && !isempty(report_mzml.matched_peaks)
                avg_mz = mean([p[1] for p in report_mzml.matched_peaks])
                mz_tolerance_mzml = avg_mz * report_mzml.optimal_ppm / 1e6
                println("Optimal PPM: $(round(report_mzml.optimal_ppm, digits=2)), Average m/z: $(round(avg_mz, digits=2))")
                println("Calculated m/z tolerance: $(round(mz_tolerance_mzml, digits=5))")
            else
                println("Could not determine optimal tolerance, using default: $mz_tolerance_mzml")
            end

            # Create a dedicated subdirectory for the output plots
            mzml_output_dir = joinpath(RESULTS_DIR, "mzml_pipeline_stages")
            mkpath(mzml_output_dir)
            
            test_full_pipeline(msi_data_mzml, MZML_SPECTRUM_ID, output_dir=mzml_output_dir, file_type_prefix="mzml", mz_tolerance=mz_tolerance_mzml)
            test_full_pipeline_on_total_spectrum(msi_data_mzml, output_dir=mzml_output_dir, file_type_prefix="mzml", mz_tolerance=mz_tolerance_mzml)
        catch e
            println("ERROR in .mzML pipeline test: $e")
            showerror(stdout, e, catch_backtrace())
        end
    else
        println("SKIPPED: File not found: $TEST_MZML_FILE")
    end

    # --- Test Case 2: Run full pipeline on a single imzML spectrum ---
    println("\n" * "="^20 * " Test Case 2: Full Pipeline on .imzML Spectrum " * "="^20)
    println("FILE: ", TEST_IMZML_FILE)
    if isfile(TEST_IMZML_FILE)
        try
            msi_data_imzml = OpenMSIData(TEST_IMZML_FILE)

            # Dynamically determine tolerance
            println("\n--- Calculating optimal tolerance for .imzML data ---")
            report_imzml = analyze_mass_accuracy(msi_data_imzml, get_common_calibration_standards(:maldi_pos))
            mz_tolerance_imzml = 0.002 # Default
            if haskey(report_imzml, :optimal_ppm) && !isnan(report_imzml.optimal_ppm) && !isempty(report_imzml.matched_peaks)
                avg_mz = mean([p[1] for p in report_imzml.matched_peaks])
                mz_tolerance_imzml = avg_mz * report_imzml.optimal_ppm / 1e6
                println("Optimal PPM: $(round(report_imzml.optimal_ppm, digits=2)), Average m/z: $(round(avg_mz, digits=2))")
                println("Calculated m/z tolerance: $(round(mz_tolerance_imzml, digits=5))")
            else
                println("Could not determine optimal tolerance, using default: $mz_tolerance_imzml")
            end
            
            # Create a dedicated subdirectory for the output plots
            imzml_output_dir = joinpath(RESULTS_DIR, "imzml_pipeline_stages")
            mkpath(imzml_output_dir)
            
            test_full_pipeline(msi_data_imzml, IMZML_COORDS, output_dir=imzml_output_dir, file_type_prefix="imzml", mz_tolerance=mz_tolerance_imzml)
            test_full_pipeline_on_total_spectrum(msi_data_imzml, output_dir=imzml_output_dir, file_type_prefix="imzml", mz_tolerance=mz_tolerance_imzml)
            # generate_qc_report(msi_data_imzml, TEST_IMZML_FILE, output_dir=imzml_output_dir)
            custom_reference_peaks = Dict(
                31.974 => "Red Phosphorus",
                432.6584 => "P13",
                464.6059 => "P15",
                526.5534 => "P17",
                650.4485 => "P21",
                774.3435 => "P25",
                898.2385 => "P29",
                950.1861 => "P31",
                1022.1336 => "P33",
                1146.0286 => "P37",
                1593.8187 => "P45",
                772.433 => "Unknown 1",
                772.5253 => "Unknown 2"
            )
            n_samples = length(msi_data_imzml.spectra_metadata)
            generate_qc_report(msi_data_imzml, TEST_IMZML_FILE, reference_peaks=custom_reference_peaks, output_dir=imzml_output_dir, sample_spectra=n_samples)
        catch e
            println("ERROR in .imzML pipeline test: $e")
            showerror(stdout, e, catch_backtrace())
        end
    else
        println("SKIPPED: File not found: $TEST_IMZML_FILE")
    end

    # --- Test Case 3: Generate QC Report for .imzML data ---
    println("\n" * "="^20 * " Test Case 3: QC Report Generation for .imzML " * "="^20)
    println("FILE: ", TEST_IMZML_FILE)
    if isfile(TEST_IMZML_FILE)
        try
            msi_data_imzml = OpenMSIData(TEST_IMZML_FILE)
            
            # Create a dedicated subdirectory for the QC report
            qc_output_dir = joinpath(RESULTS_DIR, "qc_report")
            mkpath(qc_output_dir)

            println("\n--- Generating comprehensive QC report ---")
            custom_reference_peaks = Dict(
                31.974 => "Red Phosphorus",
                432.6584 => "P13",
                464.6059 => "P15",
                526.5534 => "P17",
                650.4485 => "P21",
                774.3435 => "P25",
                898.2385 => "P29",
                950.1861 => "P31",
                1022.1336 => "P33",
                1146.0286 => "P37",
                1593.8187 => "P45",
                772.433 => "Unknown 1",
                772.5253 => "Unknown 2"
            )
            # You can control the number of spectra sampled for the QC report.
            # For the most accurate results, you can sample all spectra, but it will take longer.
            # To sample all, use: n_samples = length(msi_data_imzml.spectra_metadata)
            n_samples = length(msi_data_imzml.spectra_metadata)

            #generate_qc_report(msi_data_imzml, TEST_IMZML_FILE, output_dir=qc_output_dir)
            generate_qc_report(msi_data_imzml, TEST_IMZML_FILE, reference_peaks=custom_reference_peaks, output_dir=qc_output_dir, sample_spectra=n_samples)
            
            println("\n--- Analyzing specific reference peaks ---")
            #=
            reference_peaks = Dict(
                104.10754 => "Imidazole",
                175.11995 => "GPC fragment", 
                226.15687 => "Phosphocholine"
            )
            =#
            report = analyze_mass_accuracy(msi_data_imzml, custom_reference_peaks)
            if haskey(report, :optimal_ppm)
                println("Optimal PPM tolerance with specific peaks: $(round(report.optimal_ppm, digits=2)) ppm")
            else
                println("Could not determine optimal PPM with specific peaks.")
            end

        catch e
            println("ERROR in QC report generation test: $e")
            showerror(stdout, e, catch_backtrace())
        end
    else
        println("SKIPPED: File not found: $TEST_IMZML_FILE")
    end

    println("\nPreprocessing tests finished.")
end

# --- Execute ---
# Ensure the results directory exists
mkpath(RESULTS_DIR)
@time run_preprocessing_tests()