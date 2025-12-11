# MSI_src.jl

module MSI_src

# Export the public MSI API
export OpenMSIData, 
       GetSpectrum, 
       IterateSpectra, 
       ImportMzmlFile,
       plot_slice, 
       get_total_spectrum, 
       get_average_spectrum,
       precompute_analytics, 
       process_spectrum,
       generate_colorbar_image,
       process_image_pipeline,
       load_and_prepare_mask,
       set_global_mz_range!,
       get_global_mz_range,
       MSIData,
       _iterate_spectra_fast,
       validate_spectrum

# Export the public preprocessing & precalculations API
export run_preprocessing_analysis,
       main_precalculation,
       FeatureMatrix,
       Calibration,
       Smoothing,
       BaselineCorrection,
       Normalization,
       PeakPicking,
       PeakBinning,
       get_masked_spectrum_indices,
       detect_peaks_profile_core, 
       detect_peaks_centroid_core, 
       smooth_spectrum_core, 
       apply_baseline_correction_core, 
       apply_normalization_core, 
       bin_peaks_core, 
       PeakSelection,
       PeakAlignment,
       find_calibration_peaks_core, 
       align_peaks_lowess_core,
       MutableSpectrum,
       transform_intensity_core,
       calibrate_spectra_core

export apply_baseline_correction,
       apply_smoothing,
       apply_peak_picking,
       apply_calibration,
       apply_peak_alignment,
       apply_normalization,
       apply_peak_binning,
       apply_intensity_transformation,
       save_feature_matrix

# Include all source files directly into the main module
include("BloomFilters.jl")
include("Common.jl")
include("Platform.jl")
include("MSIData.jl")
include("ParserHelpers.jl")
include("mzML.jl")
include("imzML.jl")
include("MzmlConverter.jl")
include("Preprocessing.jl")
include("ImageProcessing.jl")
include("Precalculations.jl")
include("PreprocessingPipeline.jl")

using Setfield # For immutable struct updates


# --- Main Entry Point --- #

"""
    OpenMSIData(filepath::String; cache_size=100)

Opens a .mzML or .imzML file and prepares it for data access.

This is the main entry point for the new data access API.
"""
function OpenMSIData(filepath::String; cache_size=300, spectrum_type_map::Union{Dict{Int, Symbol}, Nothing}=nothing, use_mmap::Bool=false)
    local msi_data

    # Detect platform profile once at the entry point
    platform_profile = detect_platform_profile()

    if endswith(lowercase(filepath), ".mzml")
        msi_data = load_mzml_lazy(filepath, cache_size=cache_size, use_mmap=use_mmap)
    elseif endswith(lowercase(filepath), ".imzml")
        msi_data = load_imzml_lazy(filepath, cache_size=cache_size, use_mmap=use_mmap)
    else
        error("Unsupported file type: $filepath. Please provide a .mzML or .imzML file.")
    end

    # Apply spectrum type map if provided
    if spectrum_type_map !== nothing
        for (idx, type_symbol) in spectrum_type_map
            if 1 <= idx <= length(msi_data.spectra_metadata)
                msi_data.spectra_metadata[idx] = @set msi_data.spectra_metadata[idx].type = type_symbol
            else
                @warn "Spectrum index $idx out of bounds for spectrum_type_map. Skipping."
            end
        end
    end
    
    return msi_data
end

end # module MSI_src
