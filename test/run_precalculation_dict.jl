# test/run_precalculation_example.jl

using Printf
import Pkg

# --- Load the MSI_src Module ---
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src

# ===================================================================
# CONFIG: PLEASE FILL IN YOUR FILE PATHS HERE
# ===================================================================

# const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.mzML"
#const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/mzML/Col_1.mzML"
const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Thricoderma_etc/Imaging_interaccion_trichoderma_vs_streptomyces.mzML"

# const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.imzML"
#const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML"
const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Thricoderma_etc/Imaging_interaccion_trichoderma_vs_streptomyces.imzML"

#const MASK_ROUTE = "/home/pixel/Documents/Cinvestav_2025/JuliaMSI/public/css/masks/Stomach_DHB_uncompressed.png"
const MASK_ROUTE = ""

const reference_peaks = Dict(
    # DHB Matrix peaks (should be present)
    137.0244 => "DHB_fragment",
    155.0349 => "DHB_M+H", 
    177.0168 => "DHB_M+Na",
    
    # Common lipids in your mass range
    496.3398 => "PC_16:0_16:0",
    520.3398 => "PC_16:0_18:1", 
    760.5851 => "PC_16:0_18:1_Na",
    
    # Common contaminants
    391.2843 => "PDMS",
    413.2662 => "PDMS_Na",
    
    # Add some high mass peaks
    842.5092 => "Protein_standard",
    1045.532 => "Protein_standard",

    290.1747 => "Atropine [M+H]+", 
    304.1903 => "Scopolamine [M+H]+",
    124.0393 => "Tropine [M+H]+",
)

# ===================================================================
# HELPER FUNCTIONS FOR PRINTING
# ===================================================================

function print_header(title::String)
    println("\n" * "="^80)
    println("$(title)")
    println("="^80)
end

function check_data_range(msi_data::MSIData)
    println("\n--- Data Range Analysis ---")
    
    min_mz, max_mz = get_global_mz_range(msi_data)
    if isfinite(min_mz) && isfinite(max_mz) && min_mz < max_mz
        println("Global m/z range: [$(min_mz), $(max_mz)]")
    else
        println("Global m/z range: Not yet determined or invalid (initial: [$(min_mz), $(max_mz)])")
    end
    
    # Check a few spectra to see actual m/z values
    println("\nChecking first few spectra for actual m/z values:")
    for i in 1:min(3, length(msi_data.spectra_metadata))
        try
            mz, intensity = GetSpectrum(msi_data, i)
            if !isempty(mz)
                println("Spectrum $i: m/z range [$(minimum(mz)), $(maximum(mz))], length=$(length(mz))")
                # Print first and last few m/z values
                if length(mz) > 10
                    println("  First 5 m/z: $(mz[1:5])")
                    println("  Last 5 m/z: $(mz[end-4:end])")
                end
            end
        catch e
            println("Spectrum $i: Error - $e")
        end
    end
end

# ===================================================================
# MAIN EXAMPLE RUNNER
# ===================================================================

function run_precalculation_example()
    # --- Process mzML file ---
    print_header("Processing mzML File: $(basename(TEST_MZML_FILE))")
    if !isfile(TEST_MZML_FILE)
        println("SKIPPING: mzML file not found at $(TEST_MZML_FILE)")
    else
        try
            msi_data_mzml = @time OpenMSIData(TEST_MZML_FILE)
            check_data_range(msi_data_mzml)
            analysis_results_mzml = main_precalculation(msi_data_mzml, reference_peaks=reference_peaks)
            
            println("\n" * "*"^80)
            println("MZML PREPROCESSING ANALYSIS RESULTS")
            println("*"^80)
            for (step_name, params) in analysis_results_mzml
                println("\n--- $(uppercase(string(step_name))) Parameters ---")
                if isempty(params)
                    println("  No recommended parameters.")
                else
                    for (param_name, param_value) in params
                        println("  $param_name: $param_value")
                    end
                end
            end
            
            close(msi_data_mzml) # Close file handles
        catch e
            println("ERROR processing mzML file: $e")
            showerror(stdout, e, catch_backtrace())
        end
    end

    # --- Process imzML file ---
    print_header("Processing imzML File: $(basename(TEST_IMZML_FILE))")
    if !isfile(TEST_IMZML_FILE)
        println("SKIPPING: imzML file not found at $(TEST_IMZML_FILE)")
    else
        try
            msi_data_imzml = @time OpenMSIData(TEST_IMZML_FILE)
            check_data_range(msi_data_imzml)
            analysis_results_imzml = main_precalculation(msi_data_imzml, reference_peaks=reference_peaks, mask_path=MASK_ROUTE)
            
            println("\n" * "*"^80)
            println("IMZML PREPROCESSING ANALYSIS RESULTS")
            println("*"^80)
            for (step_name, params) in analysis_results_imzml
                println("\n--- $(uppercase(string(step_name))) Parameters ---")
                if isempty(params)
                    println("  No recommended parameters.")
                else
                    for (param_name, param_value) in params
                        println("  $param_name: $param_value")
                    end
                end
            end

            close(msi_data_imzml) # Close file handles
        catch e
            println("ERROR processing imzML file: $e")
            showerror(stdout, e, catch_backtrace())
        end
    end
end

# --- Execute ---
@time run_precalculation_example()
