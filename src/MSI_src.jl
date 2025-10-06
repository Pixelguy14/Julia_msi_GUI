module MSI_src

# Export the public API
export OpenMSIData, GetSpectrum, IterateSpectra, ImportMzmlFile, load_slices, plot_slices, plot_slice, get_total_spectrum, get_average_spectrum, LoadMzml, LoadSpectra

# Include all source files directly into the main module
include("MSIData.jl")
include("ParserHelpers.jl")
include("mzML.jl")
include("imzML.jl")
include("MzmlConverter.jl")


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
