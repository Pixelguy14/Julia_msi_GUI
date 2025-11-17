# src/PreprocessingVersioning.jl

"""
This module provides a versioning helper functions in our data preprocessing, module
inspired by the functionality of the R package MALDIquant. It includes
functions for quality control, intensity transformation, smoothing, baseline correction,
normalization, peak picking, alignment, and feature matrix generation.

This module is unfinished and may not end up in the release version.
"""

# =============================================================================
# Dependencies
# =============================================================================

using UUIDs, JLD2, CodecBase # For preprocessing versioning

# =============================================================================
# Data Structures
# =============================================================================

"""
    VersionedSpectralData

A struct to hold a snapshot of spectral data at a specific point in the
preprocessing workflow. Each transformation creates a new instance of this
struct, forming a history of all operations.

# Fields
- `version_id::String`: A unique identifier for this version of the data.
- `parent_id::String`: The identifier of the version from which this one was derived.
- `spectra::Vector`: The spectral data itself, typically a vector of `(mz, intensity)` tuples.
- `processing_step::String`: A description of the operation that created this version (e.g., "smooth", "baseline").
- `timestamp::DateTime`: The time at which this version was created.
- `parameters::Dict`: A dictionary of the parameters used in the processing step.
"""
struct VersionedSpectralData
    version_id::String
    parent_id::String
    spectra::Vector
    processing_step::String
    timestamp::DateTime
    parameters::Dict
end

"""
    MSISession

Manages the state of a preprocessing session, including the history of all
data versions and the current working version.

# Fields
- `session_id::String`: A unique identifier for the entire session.
- `original_data::VersionedSpectralData`: The initial, unprocessed spectral data.
- `history::Vector{VersionedSpectralData}`: A chronological list of all data versions created during the session.
- `current_version::VersionedSpectralData`: The version of the data currently being worked on or displayed.
- `processing_steps::Vector{String}`: A human-readable log of the processing steps applied.
"""
mutable struct MSISession
    session_id::String
    original_data::VersionedSpectralData
    history::Vector{VersionedSpectralData}
    current_version::VersionedSpectralData
    processing_steps::Vector{String}
end

# =============================================================================
# Session Management
# =============================================================================

"""
    initialize_session(raw_spectra::Vector, session_name::String) -> MSISession

Creates a new preprocessing session from raw spectral data.

# Arguments
- `raw_spectra::Vector`: A vector of `(mz, intensity)` tuples representing the initial dataset.
- `session_name::String`: A name for the session.

# Returns
- `MSISession`: A new session object initialized with the provided data.
"""
function initialize_session(raw_spectra::Vector, session_name::String)::MSISession
    original = VersionedSpectralData(
        string(uuid4()),
        "root",
        raw_spectra,
        "Original Data",
        now(),
        Dict()
    )
    return MSISession(session_name, original, [original], original, [])
end

"""
    get_processing_history(session::MSISession) -> Vector{String}

Returns a list of the processing steps applied during the session.
"""
get_processing_history(session::MSISession) = [v.processing_step for v in session.history]

"""
    undo_last_step!(session::MSISession)

Reverts the session to the state before the last processing step was applied.
This modifies the session in place.
"""
function undo_last_step!(session::MSISession)
    if length(session.history) > 1
        pop!(session.history)
        pop!(session.processing_steps)
        session.current_version = last(session.history)
    end
    return session
end

"""
    revert_to_version!(session::MSISession, version_id::String)

Reverts the session to a specific version in its history, discarding all
subsequent changes. This modifies the session in place.
"""
function revert_to_version!(session::MSISession, version_id::String)
    target_idx = findfirst(v -> v.version_id == version_id, session.history)
    if !isnothing(target_idx)
        session.history = session.history[1:target_idx]
        session.current_version = session.history[end]
        # Also truncate the descriptive processing steps
        num_steps_to_keep = max(0, target_idx - 1)
        session.processing_steps = session.processing_steps[1:num_steps_to_keep]
    end
    return session
end

# =============================================================================
# Versioned Preprocessing
# =============================================================================

"""
    apply_processing_step(spectra::Vector, step::Symbol, params::Dict) -> Vector

A dispatcher that applies a single, specified preprocessing function to a set of spectra.
This is the core function called by the versioned workflow.

# Arguments
- `spectra::Vector`: The input spectral data.
- `step::Symbol`: The symbol representing the processing step (e.g., `:smooth`, `:baseline`).
- `params::Dict`: A dictionary of parameters for the step.

# Returns
- `Vector`: The processed spectral data.
"""
function apply_processing_step(spectra::Vector, step::Symbol, params::Dict)
    processed = deepcopy(spectra)
    params_sym = Dict(Symbol(k) => v for (k,v) in params) # Ensure keys are symbols

    if step === :transform
        for i in eachindex(processed)
            mz, y = processed[i]
            processed[i] = (mz, transform_intensity(y; params_sym...))
        end
    elseif step === :smooth
        for i in eachindex(processed)
            mz, y = processed[i]
            processed[i] = (mz, smooth_spectrum(y; params_sym...))
        end
    elseif step === :baseline
        for i in eachindex(processed)
            mz, y = processed[i]
            baseline = snip_baseline(y; params_sym...)
            processed[i] = (mz, max.(0.0, y .- baseline))
        end
    elseif step === :normalize
        mode = get(params, :normalize_method, :tic)
        if mode === :pqn
            matrix = hcat([float.(p[2]) for p in processed]...)
            matrix_norm = pqn_normalize(matrix)
            for i in eachindex(processed)
                mz, _ = processed[i]
                processed[i] = (mz, view(matrix_norm, :, i))
            end
        else # :tic or :median
            norm_func = (mode === :tic) ? tic_normalize : median_normalize
            for i in eachindex(processed)
                mz, y = processed[i]
                processed[i] = (mz, norm_func(y))
            end
        end
    elseif step === :peaks
        peak_results = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, length(processed))
        for (i, (mz, y)) in enumerate(processed)
            pk_mz, pk_int = detect_peaks_profile(mz, y; params_sym...)
            peak_results[i] = (pk_mz, pk_int)
        end
        return peak_results
    elseif step === :peaks_wavelet
        peak_results = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, length(processed))
        for (i, (mz, y)) in enumerate(processed)
            pk_mz, pk_int = detect_peaks_wavelet(mz, y; params_sym...)
            peak_results[i] = (pk_mz, pk_int)
        end
        return peak_results
    elseif step === :calibrate
        return calibrate_spectra(processed; params_sym...)
    else
        @warn "Unsupported processing step: $step"
    end
    return processed
end

"""
    run_versioned_preprocessing(session::MSISession, step::Symbol, params::Dict) -> VersionedSpectralData

Creates a new version of spectral data by applying a single processing step to the
current version in the session.

# Arguments
- `session::MSISession`: The current processing session.
- `step::Symbol`: The processing step to apply.
- `params::Dict`: Parameters for the processing step.

# Returns
- `VersionedSpectralData`: A new data version with the transformation applied.
"""
function run_versioned_preprocessing(session::MSISession, step::Symbol, params::Dict)
    # Create a new version based on the current one
    parent_version = session.current_version
    new_version_id = string(uuid4())
    
    # Apply the processing step
    new_spectra = apply_processing_step(parent_version.spectra, step, params)
    
    # Create the new versioned data object
    new_version = VersionedSpectralData(
        new_version_id,
        parent_version.version_id,
        new_spectra,
        string(step),
        now(),
        params
    )
    
    return new_version
end

"""
    apply_processing_with_version!(session::MSISession, step::Symbol, ui_params::Dict)

A high-level function to apply a processing step, create a new version, and
update the session state in place.

# Arguments
- `session::MSISession`: The session to modify.
- `step::Symbol`: The processing step to apply.
- `ui_params::Dict`: A dictionary of parameters, typically from a UI.
"""
function apply_processing_with_version!(session::MSISession, step::Symbol, ui_params::Dict)
    # In a real app, you would validate and convert UI params here.
    processing_params = ui_params
    
    # Create the new version
    new_version = run_versioned_preprocessing(session, step, processing_params)
    
    # Update the session
    push!(session.history, new_version)
    session.current_version = new_version
    
    # Describe the step for the history log
    step_description = "$(string(step)) with params: " * join(["$k=$v" for (k,v) in processing_params], ", ")
    push!(session.processing_steps, step_description)
    
    return session
end

# =============================================================================
# Data Serialization & Export
# =============================================================================

"""
    save_spectral_version(data::VersionedSpectralData, filepath::String)

Saves a `VersionedSpectralData` object to a file using the JLD2 format.
"""
function save_spectral_version(data::VersionedSpectralData, filepath::String)
    JLD2.save_object(filepath, data)
end

"""
    load_spectral_version(filepath::String) -> VersionedSpectralData

Loads a `VersionedSpectralData` object from a JLD2 file.
"""
function load_spectral_version(filepath::String)::VersionedSpectralData
    return JLD2.load_object(filepath)
end

"""
    _encode_base64(data::AbstractVector{<:Real}) -> String

Helper function to convert a numeric vector into a Base64 encoded string.
"""
function _encode_base64(data::AbstractVector{T}) where T <: Real
    bytes = reinterpret(UInt8, data)
    return base64encode(bytes)
end

"""
    export_to_mzml(session::MSISession, filepath::String)

Exports the current version of the spectral data in a session to a standard
.mzML file.

# Arguments
- `session::MSISession`: The current session.
- `filepath::String`: The path for the output .mzML file.
"""
function export_to_mzml(session::MSISession, filepath::String)
    spectra_to_export = session.current_version.spectra
    
    open(filepath, "w") do f
        # XML Header
        write(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        write(f, "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" version=\"1.1\">\n")
        
        # CV List
        write(f, "  <cvList count=\"2\">\n")
        write(f, "    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.123\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n")
        write(f, "    <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"07:03:2022\" URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>\n")
        write(f, "  </cvList>\n")

        # File Description
        write(f, "  <fileDescription>\n    <fileContent>\n")
        write(f, "      <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\"/>\n")
        write(f, "    </fileContent>\n  </fileDescription>\n")

        # Run and Spectrum List
        write(f, "  <run id=\"run1\" defaultInstrumentConfigurationRef=\"instrument1\">\n")
        write(f, "    <spectrumList count=\"$(length(spectra_to_export))\" defaultDataProcessingRef=\"dp1\">\n")

        for (i, (mz, intensity)) in enumerate(spectra_to_export)
            # Ensure data is in the correct format for encoding
            mz_64 = convert(Vector{Float64}, mz)
            int_32 = convert(Vector{Float32}, intensity)

            # Base64 encode the binary data
            mz_b64 = _encode_base64(mz_64)
            int_b64 = _encode_base64(int_32)

            write(f, "      <spectrum index=\"$(i-1)\" id=\"scan=$(i)\" defaultArrayLength=\"$(length(mz))\">\n")
            write(f, "        <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>\n")
            # Assuming profile mode for processed data, could be made dynamic
            write(f, "        <cvParam cvRef=\"MS\" accession=\"MS:1000128\" name=\"profile spectrum\"/>\n")
            
            write(f, "        <binaryDataArrayList count=\"2\">\n")
            
            # m/z array
            write(f, "          <binaryDataArray encodedLength=\"$(length(mz_b64))\">\n")
            write(f, "            <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n")
            write(f, "            <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>\n")
            write(f, "            <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n")
            write(f, "            <binary>$(mz_b64)</binary>\n")
            write(f, "          </binaryDataArray>\n")

            # Intensity array
            write(f, "          <binaryDataArray encodedLength=\"$(length(int_b64))\">\n")
            write(f, "            <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>\n")
            write(f, "            <cvParam cvRef=\"MS\" accession=\"MS:1000521\" name=\"32-bit float\"/>\n")
            write(f, "            <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>\n")
            write(f, "            <binary>$(int_b64)</binary>\n")
            write(f, "          </binaryDataArray>\n")
            
            write(f, "        </binaryDataArrayList>\n")
            write(f, "      </spectrum>\n")
        end

        write(f, "    </spectrumList>\n")
        write(f, "  </run>\n")
        
        # Dummy instrument and data processing info
        write(f, "  <instrumentConfigurationList count=\"1\">\n")
        write(f, "    <instrumentConfiguration id=\"instrument1\"/>\n")
        write(f, "  </instrumentConfigurationList>\n")
        write(f, "  <dataProcessingList count=\"1\">\n")
        write(f, "    <dataProcessing id=\"dp1\"/>\n")
        write(f, "  </dataProcessingList>\n")

        write(f, "</mzML>\n")
    end
    println("Successfully exported current data to $filepath")
end


# =============================================================================
# Stubs for UI Integration
# =============================================================================

"""
    validate_ui_parameters(step::Symbol, ui_params::Dict) -> Tuple{Bool, String}

(Stub) Validates parameters from a UI before they are used in a processing step.
"""
function validate_ui_parameters(step::Symbol, ui_params::Dict)::Tuple{Bool, String}
    # In a real implementation, this would check types, ranges, etc.
    # For example, for :smooth, ensure 'sg_window' is an odd integer.
    println("Validating parameters for step: $step")
    return (true, "Parameters are valid.")
end

"""
    run_processing_with_progress(session, step, params, progress_callback)

(Stub) A wrapper for running a processing step that includes a progress reporting callback.
"""
function run_processing_with_progress(session, step, params, progress_callback)
    progress_callback(0.0, "Starting $step...")
    
    # This is a simplified example. Real implementation would need to
    # hook into the loops inside apply_processing_step.
    new_version = run_versioned_preprocessing(session, step, params)
    
    progress_callback(1.0, "Finished $step.")
    return new_version
end

"""
    safe_processing_application!(session::MSISession, step::Symbol, params::Dict)

(Stub) A safe wrapper to apply a processing step that includes error handling and rollback.
"""
function safe_processing_application!(session::MSISession, step::Symbol, params::Dict)
    num_history = length(session.history)
    try
        println("Safely applying step: $step")
        apply_processing_with_version!(session, step, params)
    catch e
        @error "Processing step '$step' failed!" exception=(e, catch_backtrace())
        # Roll back to the previous state
        while length(session.history) > num_history
            undo_last_step!(session)
        end
        println("Session has been rolled back to the previous state.")
        return session # Return the rolled-back session
    end
    return session # Return the updated session
end

