# MSI_src.jl

module MSI_src

# Export the public MSI API
export OpenMSIData, 
       GetSpectrum, 
       IterateSpectra, 
       ImportMzmlFile, 
       load_slices, 
       plot_slices, 
       plot_slice, 
       get_total_spectrum, 
       get_average_spectrum, 
       LoadMzml, 
       precompute_analytics, 
       process_spectrum,
       generate_colorbar_image,
       process_image_pipeline,
       load_and_prepare_mask

# Export the public Preprocessing API
export FeatureMatrix,
       run_preprocessing_pipeline,
       qc_is_empty,
       qc_is_regular,
       transform_intensity,
       smooth_spectrum,
       snip_baseline,
       tic_normalize,
       pqn_normalize,
       detect_peaks_profile,
       align_peaks_lowess,
       bin_peaks,
       plot_stage_spectrum,
       calculate_ppm_error, 
       calculate_resolution_fwhm, 
       analyze_mass_accuracy,
       generate_qc_report, 
       get_common_calibration_standards

# Include all source files directly into the main module
include("MSIData.jl")
include("ParserHelpers.jl")
include("mzML.jl")
include("imzML.jl")
include("MzmlConverter.jl")
include("Preprocessing.jl")
include("ImageProcessing.jl")


# --- Main Entry Point --- #

"""
    OpenMSIData(filepath::String; cache_size=100)

Opens a .mzML or .imzML file and prepares it for data access.

This is the main entry point for the new data access API.
"""
function OpenMSIData(filepath::String; cache_size=300)
    if endswith(lowercase(filepath), ".mzml")
        return load_mzml_lazy(filepath, cache_size=cache_size)
    elseif endswith(lowercase(filepath), ".imzml")
        return load_imzml_lazy(filepath, cache_size=cache_size)
    else
        error("Unsupported file type: $filepath. Please provide a .mzML or .imzML file.")
    end
end

end # module MSI_src
