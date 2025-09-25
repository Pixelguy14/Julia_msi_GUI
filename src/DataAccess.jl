# src/DataAccess.jl

"""
    load_spectra(fileName::String)

Eagerly loads all spectra from a .mzML or .imzML file into a 2xN matrix.
This function follows the eager-loading pattern of the original `LoadMzml`
for backward compatibility and for analyses that require all data in memory.

Each column in the returned matrix represents a single spectrum:
- Row 1: m/z array
- Row 2: Intensity array

For .imzML files, the (x,y) spatial coordinates are discarded.

# Arguments
* `fileName`: Full path to the .mzML or .imzML file.

# Returns
- An `Array{Any, 2}` containing all spectra.
"""
function load_spectra(fileName::String)
    if endswith(lowercase(fileName), ".mzml")
        # Detected .mzML file. Using the existing eager loader.
        return LoadMzml(fileName)

    elseif endswith(lowercase(fileName), ".imzml")
        # Detected .imzML file. Eagerly loading all spectra from .ibd.
        
        # Lazily load metadata to get offsets and counts
        imzml_data = load_imzml(fileName)
        
        try
            num_spectra = length(imzml_data.spectra_metadata)
            if num_spectra == 0
                return Array{Any}(undef, (2, 0))
            end

            spectra_matrix = Array{Any}(undef, (2, num_spectra))
            
            hIbd = imzml_data.ibd_handle
            mz_format = imzml_data.mz_format
            int_format = imzml_data.intensity_format

            # For performance, pre-allocate one large buffer for reading.
            # This avoids re-allocating memory for each spectrum inside the loop.
            max_points = maximum(meta.mz_count for meta in imzml_data.spectra_metadata)
            mz_array_buffer = Array{mz_format}(undef, max_points)
            intensity_array_buffer = Array{int_format}(undef, max_points)

            for i in 1:num_spectra
                meta = imzml_data.spectra_metadata[i]
                
                # Create a view into the buffer with the correct size for this spectrum
                current_mz_array = view(mz_array_buffer, 1:meta.mz_count)
                current_int_array = view(intensity_array_buffer, 1:meta.intensity_count)

                # Read binary data directly into the sized views
                seek(hIbd, meta.mz_offset)
                read!(hIbd, current_mz_array)
                
                seek(hIbd, meta.intensity_offset)
                read!(hIbd, current_int_array)

                # Store a copy in the final matrix. A copy is necessary because
                # the buffer is overwritten in the next iteration.
                spectra_matrix[1, i] = copy(current_mz_array)
                spectra_matrix[2, i] = copy(current_int_array)
            end
            
            return spectra_matrix
        finally
            # Ensure the .ibd file handle is closed, as this is an eager load.
            if isopen(imzml_data.ibd_handle)
                close(imzml_data.ibd_handle)
            end
        end

    else
        error("Unsupported file type. Please provide a .mzML or .imzML file.")
    end
end
