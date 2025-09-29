# test/run_tests.jl

# ===================================================================
# Test Environment for JuliaMSI Package
# ===================================================================
# This script validates the core functionality of the data processing
# workflows, including loading, converting, and visualizing mass
# spectrometry data.
#
# Instructions:
# 1. Fill in the placeholder paths in the "CONFIG" section below.
# 2. Run the script from the project's root directory:
#    julia test/run_tests.jl
# 3. Check the `test/results/` folder for the output images.
# ===================================================================

using Printf
using CairoMakie
import Pkg

# --- Load the MSI_src Module ---
# Activate the project environment at the parent directory of this test script
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src


# ===================================================================
# CONFIG: PLEASE FILL IN YOUR FILE PATHS HERE
# ===================================================================

# --- Test Case 1: Standard .mzML file ---
# A regular, non-imaging mzML file.
# const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/mzML/T9_A1.mzML"
# const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.mzML"
const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.mzML"
const SPECTRUM_TO_PLOT = 1 # Which spectrum to plot from the file

# --- Test Case 2: .mzML + Sync File for Conversion ---
# The special .mzML file with one spectrum per pixel.
const CONVERSION_SOURCE_MZML = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.mzML"
# The corresponding synchronization text file.
const CONVERSION_SYNC_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.txt"

# const CONVERSION_SOURCE_MZML = "/home/pixel/Documents/Cinvestav_2025/Analisis/Imaging_paper_spray/Imaging_paper_spray.mzML"
# const CONVERSION_SYNC_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Imaging_paper_spray/Imaging_paper_spray.txt"

# The desired output path for the new .imzML file.
const CONVERSION_TARGET_IMZML = "test/results/converted_mzml.imzML"

# --- Test Case 3: Standard .imzML file ---
# An existing imzML file (can be the one generated from Case 2).
const TEST_IMZML_FILE = CONVERSION_TARGET_IMZML # The output from case 2
# const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML"
# const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Imaging_paper_spray/Imaging_paper_spray.imzML"
# The m/z value to use for creating an image slice.
const MZ_VALUE_FOR_SLICE = 309.06 # BF
# const MZ_VALUE_FOR_SLICE = 896.0 # HR2MSI
# const MZ_VALUE_FOR_SLICE = 76.03 # I PS
# const MZ_TOLERANCE = 0.1
const MZ_TOLERANCE = 1

# Coordinates to plot a specific spectrum from imzML
const COORDS_TO_PLOT = (50, 50) # Example coordinates (X, Y)

# --- Output Directory ---
const RESULTS_DIR = "test/results"

test1 = true
test2 = true 
test3 = true

# ===================================================================
# DATA VALIDATION UTILITY
# ===================================================================

"""
    validate_msi_data(filepath::String)

Performs a series of checks on a .mzML or .imzML file using the MSIData API.
"""
function validate_msi_data(filepath::String)
    println("" * "-"^10 * " Running validation for $filepath " * "-"^10)
    if !isfile(filepath)
        println("SKIPPED VALIDATION: File not found: $filepath")
        return false
    end

    try
        # 1. Basic structure validation
        println("Opening file with OpenMSIData...")
        msi_data = @time OpenMSIData(filepath)
        
        # 2. Compare spectrum counts
        num_spectra = length(msi_data.spectra_metadata)
        println("Found $num_spectra spectra")
        @assert num_spectra > 0 "No spectra found in file."

        # 3. Test random access
        println("Testing random access to spectra...")
        test_indices = unique([1, max(1, num_spectra ÷ 2), num_spectra])
        println("Testing indices: $test_indices")
        for idx in test_indices
            print("Fetching spectrum #$idx... ")
            mz, intensity = @time GetSpectrum(msi_data, idx)
            @assert length(mz) == length(intensity) "Spectrum $idx: mz/intensity length mismatch. Got $(length(mz)) mz values and $(length(intensity)) intensity values."
            println("OK, $(length(mz)) points.")
        end
        
        # 4. Test iteration
        println("Testing iteration over all spectra...")
        count = 0
        iter_time = @elapsed for (idx, (mz, intensity)) in IterateSpectra(msi_data)
            count += 1
            # Basic data validation
            @assert all(isfinite, mz) "Non-finite mz values in spectrum $idx"
            @assert all(>=(0), intensity) "Negative intensities in spectrum $idx"
        end
        println("Iterated over $count spectra in $iter_time seconds.")
        @assert count == num_spectra "Iteration count mismatch: expected $num_spectra, got $count."
        
        println("VALIDATION SUCCESSFUL for $filepath")
        return true

    catch e
        println("VALIDATION FAILED for $filepath.")
        showerror(stdout, e, catch_backtrace())
        println()
        return false
    end
end


# ===================================================================
# TEST RUNNER
# ===================================================================

function run_test()
    println("Starting MSI_src Test Suite...")
    
    # --- Test Case 1: Process a standard .mzML file ---
    println("" * "="^20 * " Test Case 1: Processing .mzML " * "="^20)
    if test1 == true
        # Run new, stronger validation
        validate_msi_data(TEST_MZML_FILE)

        # Also run original plotting test to ensure visualization still works
        if isfile(TEST_MZML_FILE)
            try
                # Get the msi data from the mzml
                println("Plotting a sample spectrum from $TEST_MZML_FILE...")
                msi_data = OpenMSIData(TEST_MZML_FILE)
                mz, intensity = GetSpectrum(msi_data, SPECTRUM_TO_PLOT)
                
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Intensity", title="Spectrum #$SPECTRUM_TO_PLOT from $(basename(TEST_MZML_FILE))")
                lines!(ax, mz, intensity)
                
                output_path = joinpath(RESULTS_DIR, "test_mzml_spectrum.png")
                save(output_path, fig)
                println("SUCCESS: Spectrum plot saved to $output_path")
                # Get the summed spectrum data
                mz, intensity = get_total_spectrum(msi_data)

                # Plot the data
                println("Plotting total spectrum...")
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Total Intensity", title="Total Spectrum from $(basename(TEST_MZML_FILE))")
                lines!(ax, mz, intensity)

                # Saving the output
                output_path = joinpath(RESULTS_DIR, "test_mzml_total_spectrum.png")
                save(output_path, fig)
                println("SUCCESS: Total spectrum plot saved to $output_path")

                # Get the averaged spectrum data
                mz, intensity = get_average_spectrum(msi_data)

                # Plot the data
                println("Plotting averaged spectrum...")
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Average Intensity", title="Average Spectrum from $(basename(TEST_MZML_FILE))")
                lines!(ax, mz, intensity)

                # Saving the output
                output_path = joinpath(RESULTS_DIR, "test_mzml_average_spectrum.png")
                save(output_path, fig)
                println("SUCCESS: Total spectrum plot saved to $output_path")
            catch e
                println("ERROR during plotting in Test Case 1: $e")
            end
        end
    else
        println("SKIPPED Test Case 1.")
    end

    # --- Test Case 2: Convert .mzML + .txt to .imzML ---
    println("" * "="^20 * " Test Case 2: Converting to .imzML " * "="^20)
    if isfile(CONVERSION_SOURCE_MZML) && isfile(CONVERSION_SYNC_FILE) && test2 == true
        try
            println("Running conversion process (with profiling)...")
            success = @time ImportMzmlFile(CONVERSION_SOURCE_MZML, CONVERSION_SYNC_FILE, CONVERSION_TARGET_IMZML)
            if success
                println("SUCCESS: Conversion process completed.")
                # Validate the newly created imzML file
                validate_msi_data(CONVERSION_TARGET_IMZML)
            else
                println("FAILURE: Conversion process failed.")
            end
        catch e
            println("ERROR in Test Case 2: $e")
        end
    else
        println("SKIPPED: Files not found for Test Case 2.")
        println("  - mzML: $CONVERSION_SOURCE_MZML")
        println("  - Sync: $CONVERSION_SYNC_FILE")
    end

    # --- Test Case 3: Process an existing .imzML file ---
    println("
" * "="^20 * " Test Case 3: Processing .imzML " * "="^20)
    if test3 == true
        # Run new, stronger validation for imzML
        validate_msi_data(TEST_IMZML_FILE)

        # Also run tests for plotting spectrum and image slice
        if isfile(TEST_IMZML_FILE)
            # Add spectrum plotting for imzML to match Test Case 1
            try
                # Get the msi data from the imzml
                println("Plotting a sample spectrum from $TEST_IMZML_FILE...")
                msi_data = OpenMSIData(TEST_IMZML_FILE)
                x_coord, y_coord = COORDS_TO_PLOT
                
                # Get the x y coordinate spectrum data
                mz, intensity = GetSpectrum(msi_data, x_coord, y_coord)
                
                # Plot the data
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Intensity", title="Spectrum at ($x_coord, $y_coord) from $(basename(TEST_IMZML_FILE))")
                lines!(ax, mz, intensity)
                
                # Saving the output
                output_path = joinpath(RESULTS_DIR, "test_imzml_spectrum.png")
                save(output_path, fig)
                println("SUCCESS: Spectrum plot saved to $output_path")

                # Get the summed spectrum data
                mz, intensity = get_total_spectrum(msi_data)

                # Plot the data
                println("Plotting total spectrum...")
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Total Intensity", title="Total Spectrum from $(basename(TEST_IMZML_FILE))")
                lines!(ax, mz, intensity)

                # Saving the output
                output_path = joinpath(RESULTS_DIR, "test_imzml_total_spectrum.png")
                save(output_path, fig)
                println("SUCCESS: Total spectrum plot saved to $output_path")

                # Get the averaged spectrum data
                mz, intensity = get_average_spectrum(msi_data)

                # Plot the data
                println("Plotting averaged spectrum...")
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xlabel="m/z", ylabel="Average Intensity", title="Average Spectrum from $(basename(TEST_IMZML_FILE))")
                lines!(ax, mz, intensity)

                # Saving the output
                output_path = joinpath(RESULTS_DIR, "test_imzml_average_spectrum.png")
                save(output_path, fig)
                println("SUCCESS: Total spectrum plot saved to $output_path")
            catch e
                println("ERROR during spectrum plotting in Test Case 3: $e")
            end

            # Test the plot_slice function
            try
                println("Testing plot_slice function on $TEST_IMZML_FILE...")
                msi_data = OpenMSIData(TEST_IMZML_FILE)
                plot_slice(msi_data, MZ_VALUE_FOR_SLICE, MZ_TOLERANCE, RESULTS_DIR, stage_name="test_imzml_single_slice")
                # The success message is now inside plot_slice
            catch e
                println("ERROR during plot_slice test in Test Case 3: $e")
            end
        end
    else
        println("SKIPPED Test Case 3.")
    end
    
    println("Tests for all 3 cases is finished.")

end

# --- Execute ---
# Ensure the results directory exists
mkpath(RESULTS_DIR)
@time run_test()