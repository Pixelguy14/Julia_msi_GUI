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
# const TEST_MZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/mzML/Col_1.mzML"
const TEST_MZML_FILE = ""
# const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/CE4_BF_R1/CE4_BF_R1.imzML"
#const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML"
const TEST_IMZML_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Thricoderma_etc/Imaging_interaccion_trichoderma_vs_streptomyces.imzML"
#const MASK_ROUTE = "/home/pixel/Documents/Cinvestav_2025/JuliaMSI/public/css/masks/Stomach_DHB_uncompressed.png"
const MASK_ROUTE = ""

#=
const reference_peaks = Dict(
    # Common ESI positive mode reference compounds
    121.0509 => "Purine",
    149.0233 => "HP-921",
    322.0481 => "Hexakis(1H,1H,3H-tetrafluoropropoxy)phosphazine",
    622.0290 => "Hexakis(2,2-difluoroethoxy)phosphazine",
    
    # Atropine and related compounds in positive mode
    290.1747 => "Atropine [M+H]+", 
    304.1903 => "Scopolamine [M+H]+",
    124.0393 => "Tropine [M+H]+",
    
    # Common contaminants and lock masses
    391.2843 => "Polydimethylcyclosiloxane [M+H]+",
    413.2662 => "Polydimethylcyclosiloxane [M+Na]+",
    429.2402 => "Polydimethylcyclosiloxane [M+K]+"
)
=#

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
)

# ===================================================================
# HELPER FUNCTIONS FOR PRINTING
# ===================================================================

function print_header(title::String)
    println("\n" * "="^80)
    println("$(title)")
    println("="^80)
end

function print_subheader(title::String)
    println("\n" * "-"^80)
    println("$(title)")
    println("-"^80)
end

function print_param(key, value)
    if value === nothing || (isa(value, Number) && isnan(value))
        println("  " * "" * rpad(key, 30) * ": unable to get, user needs to input manually")
    elseif isa(value, Symbol) && startswith(string(key), "method") # Heuristic for method selection
        println("  " * "" * rpad(key, 30) * ": automatically detected this as the most optimal method: $(value)")
    else
        println("  " * "" * rpad(key, 30) * ": $(value)")
    end
end

function print_recommendations(recommendations::Dict)
    for (step, params) in recommendations
        print_subheader("Recommendations for $(String(step))")
        for (key, value) in params
            print_param(key, value)
        end
    end
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
            analysis_results_mzml = run_preprocessing_analysis(msi_data_mzml, reference_peaks=reference_peaks)
            
            println("\n" * "*"^80)
            println("MZML PREPROCESSING ANALYSIS RESULTS")
            println("*"^80)
            
            for (phase, results) in analysis_results_mzml
                if phase == :recommendations
                    print_subheader("Generated Preprocessing Recommendations")
                    print_recommendations(results)
                else
                    print_subheader("Phase: $(String(phase))")
                    if isa(results, Dict)
                        for (key, value) in results
                            print_param(key, value)
                        end
                    elseif isa(results, NamedTuple)
                        for field in fieldnames(typeof(results))
                            print_param(field, getfield(results, field))
                        end
                    else
                        println("  $(String(phase)) results: $(results)")
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
            analysis_results_imzml = run_preprocessing_analysis(msi_data_imzml, reference_peaks=reference_peaks, mask_path=MASK_ROUTE)
            
            println("\n" * "*"^80)
            println("IMZML PREPROCESSING ANALYSIS RESULTS")
            println("*"^80)

            for (phase, results) in analysis_results_imzml
                if phase == :recommendations
                    print_subheader("Generated Preprocessing Recommendations")
                    print_recommendations(results)
                else
                    print_subheader("Phase: $(String(phase))")
                    if isa(results, Dict)
                        for (key, value) in results
                            print_param(key, value)
                        end
                    elseif isa(results, NamedTuple)
                        for field in fieldnames(typeof(results))
                            print_param(field, getfield(results, field))
                        end
                    else
                        println("  $(String(phase)) results: $(results)")
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

# =============================================================================
# Example Usage
# =============================================================================

#=
"""
    example_preanalysis_workflow(msi_data::MSIData)

Demonstrates how to run the pre-analysis pipeline.
"""
function example_preanalysis_workflow(msi_data::MSIData)
    # Load your MSIData object (this would come from your actual data loading)
    # msi_data = load_imzml_dataset("path/to/your/data.imzML") # Assuming msi_data is passed
    
    # Define reference peaks for mass accuracy analysis
    reference_peaks = Dict(
        89.04767 => "Alanin",
        147.07642 => "Lysin",
        189.12392 => "Unknown",
        524.26496 => "PC(34:1) [M+H]+"
    )
    
    # Define region masks if you have spatial annotations
    region_masks = Dict{Symbol, BitMatrix}()
    # region_masks[:tumor] = load_mask("tumor_mask.png") # Not defined, comment out
    # region_masks[:stroma] = load_mask("stroma_mask.png") # Not defined, comment out
    
    println("Starting comprehensive pre-analysis...")
    
    # Run the complete pre-analysis pipeline
    analysis_results = run_preprocessing_analysis(
        msi_data,  # Your MSIData object
        reference_peaks=reference_peaks,
        region_masks=region_masks,
        sample_size=200  # Adjust based on dataset size
    )
    
    # Access the recommendations
    recommendations = analysis_results[:recommendations]
    
    println("\n" * "="^60)
    println("PREPROCESSING RECOMMENDATIONS")
    println("="^60)
    
    for (step, params) in recommendations
        println("\n$step:")
        for (key, value) in params
            println("  - $key: $value")
        end
    end
    
    return analysis_results
end

# You can also run individual analysis steps:
function run_targeted_analysis(msi_data::MSIData)
    # Just analyze signal quality and peak characteristics
    signal_analysis = analyze_signal_quality(msi_data, sample_size=100)
    peak_analysis = analyze_peak_characteristics(msi_data, sample_size=50)
    
    return (signal_analysis, peak_analysis)
end

=#
