# src/PreprocessingPipeline.jl

using Base.Threads # For multithreading
using Printf # For @sprintf
using Interpolations # For linear_interpolation
using DataFrames # For saving feature matrix
using CSV # For saving feature matrix

# This file provides a set of functions that apply preprocessing steps to a vector of
# `MutableSpectrum` objects. Each function takes the vector of spectra and a dictionary
# of parameters, modifying the spectra in-place where appropriate. This mirrors the
# logic from `test/run_preprocessing.jl` but is intended for use in the main application.

# ===================================================================
# PREPROCESSING PIPELINE FUNCTIONS (IN-PLACE)
# ===================================================================

"""
    apply_baseline_correction(spectra::Vector{MutableSpectrum}, params::Dict)

Applies baseline correction to the intensity data of each spectrum. This function
modifies the `.intensity` field of each `MutableSpectrum` object in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The algorithm to use. Supports `:snip`, `:convex_hull`, `:median`. Defaults to `:snip`.
- `:iterations` (Int): The number of iterations for the SNIP algorithm. Defaults to 100.
- `:window` (Int): The window size for the Median algorithm. Defaults to 20.
"""
function apply_baseline_correction(spectra::Vector{MutableSpectrum}, params::Dict)
    method = get(params, :method, :snip)
    iterations = get(params, :iterations, 100)
    window = get(params, :window, 20)

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            baseline = apply_baseline_correction_core(s.intensity; method=method, iterations=iterations, window=window)
            s.intensity = max.(0.0, s.intensity .- baseline)
        end
    end
end

"""
    apply_intensity_transformation(spectra::Vector{MutableSpectrum}, params::Dict)

Applies an intensity transformation to the intensity data of each spectrum. This function
modifies the `.intensity` field of each `MutableSpectrum` object in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The transformation to apply. Supports `:sqrt`, `:log`, `:log2`, `:log10`, `:log1p`. Defaults to `:sqrt`.
"""
function apply_intensity_transformation(spectra::Vector{MutableSpectrum}, params::Dict)
    method = get(params, :method, :sqrt)

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            s.intensity = transform_intensity_core(s.intensity; method=method)
        end
    end
end

"""
    apply_smoothing(spectra::Vector{MutableSpectrum}, params::Dict)

Applies a smoothing filter to the intensity data of each spectrum. This function
modifies the `.intensity` field of each `MutableSpectrum` object in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The smoothing algorithm. Supports `:savitzky_golay`, `:moving_average`. Defaults to `:savitzky_golay`.
- `:window` (Int): The size of the smoothing window. Defaults to 9.
- `:order` (Int): The polynomial order for the Savitzky-Golay filter. Defaults to 2.
"""
function apply_smoothing(spectra::Vector{MutableSpectrum}, params::Dict)
    method = get(params, :method, :savitzky_golay)
    window = get(params, :window, 9)
    order = get(params, :order, 2)

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            smoothed_intensity = max.(0.0, smooth_spectrum_core(s.intensity; method=method, window=window, order=order))
            s.intensity = smoothed_intensity
        end
    end
end

"""
    apply_peak_picking(spectra::Vector{MutableSpectrum}, params::Dict)

Detects peaks in each spectrum and stores them in the `.peaks` field of each
`MutableSpectrum` object, modifying it in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The peak detection algorithm. Supports `:profile`, `:wavelet`, `:centroid`. Defaults to `:profile`.
- `:snr_threshold` (Float64): Signal-to-Noise Ratio threshold.
- `:half_window` (Int): Half-window size for local maxima detection.
- `:min_peak_prominence` (Float64): Minimum required prominence for a peak.
- `:merge_peaks_tolerance` (Float64): m/z tolerance to merge adjacent peaks.
"""
function apply_peak_picking(spectra::Vector{MutableSpectrum}, params::Dict)
    method = get(params, :method, :profile)
    snr_threshold = get(params, :snr_threshold, 3.0)
    half_window = get(params, :half_window, 10)
    min_peak_prominence = get(params, :min_peak_prominence, 0.1)
    merge_peaks_tolerance = get(params, :merge_peaks_tolerance, 0.002)

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            if method == :profile
                s.peaks = detect_peaks_profile_core(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window, min_peak_prominence=min_peak_prominence, merge_peaks_tolerance=merge_peaks_tolerance)
            elseif method == :wavelet
                s.peaks = detect_peaks_wavelet_core(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window)
            elseif method == :centroid
                s.peaks = detect_peaks_centroid_core(s.mz, s.intensity; snr_threshold=snr_threshold)
            else
                s.peaks = detect_peaks_profile_core(s.mz, s.intensity; snr_threshold=snr_threshold, half_window=half_window)
            end
        else
            s.peaks = []
        end
    end
end

"""
    apply_peak_selection(spectra::Vector{MutableSpectrum}, params::Dict)

Filters peaks within each spectrum based on quality criteria. This step removes peaks
that do not meet the specified thresholds for signal-to-noise ratio (SNR),
full width at half maximum (FWHM), and peak shape.

This function modifies the `.peaks` field of each `MutableSpectrum` object in the `spectra` vector in-place.

# Parameters from `params` Dict:
- `:min_snr` (Float64): The minimum Signal-to-Noise Ratio required for a peak to be kept.
- `:min_fwhm_ppm` (Float64): The minimum FWHM (in ppm) for a peak.
- `:max_fwhm_ppm` (Float64): The maximum FWHM (in ppm) for a peak.
- `:min_shape_r2` (Float64): The minimum R² value from a Gaussian fit, measuring peak shape quality.
"""
function apply_peak_selection(spectra::Vector{MutableSpectrum}, params::Dict)
    min_snr = get(params, :min_snr, 0.0)
    min_fwhm = get(params, :min_fwhm_ppm, 0.0)
    max_fwhm = get(params, :max_fwhm_ppm, Inf)
    min_r2 = get(params, :min_shape_r2, 0.0)

    # Handle `nothing` from params, which can happen if pre-calculation fails.
    min_snr = isnothing(min_snr) ? 0.0 : min_snr
    min_fwhm = isnothing(min_fwhm) ? 0.0 : min_fwhm
    max_fwhm = isnothing(max_fwhm) ? Inf : max_fwhm
    min_r2 = isnothing(min_r2) ? 0.0 : min_r2

    Threads.@threads for s in spectra
        if !isempty(s.peaks)
            filter!(p -> 
                p.snr >= min_snr &&
                (min_fwhm <= p.fwhm <= max_fwhm) &&
                p.shape_r2 >= min_r2,
                s.peaks
            )
        end
    end
end

"""
    apply_calibration(spectra::Vector{MutableSpectrum}, params::Dict, reference_peaks::Dict)

Performs mass calibration on each spectrum using a list of internal standards.
This function modifies the `.mz` axis of each `MutableSpectrum` object in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The calibration method. Only `:internal_standards` is currently meaningful.
- `:ppm_tolerance` (Float64): The tolerance in PPM for matching detected peaks to reference masses.
- `:fit_order` (Int): The polynomial order for the calibration fit (not yet used in this implementation, defaults to linear).
"""
function apply_calibration(spectra::Vector{MutableSpectrum}, params::Dict, reference_peaks::Dict)
    method = get(params, :method, :none)
    ppm_tolerance = get(params, :ppm_tolerance, 20.0)

    if method == :none || isempty(reference_peaks)
        return
    end

    reference_masses = collect(keys(reference_peaks))

    Threads.@threads for i in 1:length(spectra)
        s = spectra[i]
        if validate_spectrum(s.mz, s.intensity)
            matched_peaks = find_calibration_peaks_core(s.mz, s.intensity, reference_masses; ppm_tolerance=ppm_tolerance)
            if length(matched_peaks) >= 2
                measured = sort(collect(values(matched_peaks)))
                theoretical = sort(collect(keys(matched_peaks)))
                itp = linear_interpolation(measured, theoretical, extrapolation_bc=Line())
                s.mz = itp(s.mz) # Modify mz-axis in-place
            else
                @warn "Spectrum $(s.id): insufficient reference peaks ($(length(matched_peaks)) found), skipping calibration."
            end
        end
    end
end

"""
    apply_peak_alignment(spectra::Vector{MutableSpectrum}, params::Dict)

Aligns the m/z axis of all spectra to a chosen reference spectrum. This function
modifies both the `.mz` axis and the m/z values within the `.peaks` field of each
`MutableSpectrum` object in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The alignment algorithm. Supports `:lowess`, `:linear`, `:ransac`.
- `:tolerance` (Float64): The tolerance for matching peaks between spectra.
- `:tolerance_unit` (Symbol): The unit for tolerance, `:mz` or `:ppm`.
"""
function apply_peak_alignment(spectra::Vector{MutableSpectrum}, params::Dict)
    method = get(params, :method, :none)
    tolerance = get(params, :tolerance, 0.002)
    tolerance_unit = get(params, :tolerance_unit, :mz)

    if method == :none
        return
    end

    ref_find_idx = findfirst(s -> !isempty(s.peaks), spectra)
    if ref_find_idx === nothing
        @warn "Insufficient spectra with peaks for alignment. Skipping."
        return
    end

    ref_spectrum = spectra[ref_find_idx]
    ref_peaks_mz = [p.mz for p in ref_spectrum.peaks]

    Threads.@threads for s in spectra
        if s.id == ref_spectrum.id || isempty(s.peaks)
            continue
        end
        
        current_peaks_mz = [p.mz for p in s.peaks]
        alignment_func = align_peaks_lowess_core(ref_peaks_mz, current_peaks_mz; method=method, tolerance=tolerance, tolerance_unit=tolerance_unit)
        
        s.mz = alignment_func.(s.mz) # Update m/z axis
        
        # Update peak m/z values
        for i in 1:length(s.peaks)
            old_peak = s.peaks[i]
            aligned_peak_mz = alignment_func(old_peak.mz)
            s.peaks[i] = (mz=aligned_peak_mz, intensity=old_peak.intensity, fwhm=old_peak.fwhm, shape_r2=old_peak.shape_r2, snr=old_peak.snr, prominence=old_peak.prominence)
        end
    end
end

"""
    apply_normalization(spectra::Vector{MutableSpectrum}, params::Dict)

Applies intensity normalization to each spectrum. This function modifies the
`.intensity` field of each `MutableSpectrum` object in-place.

# Parameters from `params` Dict:
- `:method` (Symbol): The normalization method. Supports `:tic`, `:median`, `:rms`, `:none`.
"""
function apply_normalization(spectra::Vector{MutableSpectrum}, params::Dict)
    method = get(params, :method, :tic)

    Threads.@threads for s in spectra
        if validate_spectrum(s.mz, s.intensity)
            s.intensity = apply_normalization_core(s.intensity; method=method)
        end
    end
end
function apply_peak_binning(spectra::Vector{MutableSpectrum}, params::Dict)
    tolerance = get(params, :tolerance, 20.0)
    tolerance_unit = get(params, :tolerance_unit, :ppm)
    min_peak_per_bin = get(params, :min_peak_per_bin, 3)

    if isempty(spectra) || all(s -> isempty(s.peaks), spectra)
        @warn "No peaks found for binning. Returning empty feature matrix."
        return nothing, nothing
    end
        
    all_peaks = Vector{Tuple{Float64, Float64}}()
    for s in spectra
        for p in s.peaks
            push!(all_peaks, (p.mz, p.intensity))
        end
    end
    
    if isempty(all_peaks)
        @warn "No peaks collected for binning."
        return nothing, nothing
    end
    
    sort!(all_peaks, by=x->x[1])
    
    bin_centers = Float64[]
    bin_intensities = Float64[]
    
    i = 1
    while i <= length(all_peaks)
        current_bin_start = i
        current_peak = all_peaks[i]
        
        j = i + 1
        while j <= length(all_peaks)
            next_peak = all_peaks[j]
            tol = (tolerance_unit == :ppm) ? (current_peak[1] * tolerance / 1e6) : tolerance
            
            if (next_peak[1] - current_peak[1]) <= tol
                j += 1
            else
                break
            end
        end
        
        current_bin_end = j - 1
        bin_size = current_bin_end - current_bin_start + 1
        
        if bin_size >= min_peak_per_bin
            bin_peaks = all_peaks[current_bin_start:current_bin_end]
            
            mz_sum = sum(p[1] for p in bin_peaks)
            intensity_sum = sum(p[2] for p in bin_peaks)
            
            mz_center = mz_sum / bin_size
            avg_intensity = intensity_sum / bin_size
            
            push!(bin_centers, mz_center)
            push!(bin_intensities, avg_intensity)
        end
        
        i = j
    end
    
    if !isempty(bin_centers)
        n_bins = length(bin_centers)
        feature_matrix = Matrix{Float64}(undef, 2, n_bins)
        
        for i in 1:n_bins
            feature_matrix[1, i] = bin_centers[i]
            feature_matrix[2, i] = bin_intensities[i]
        end
        
        bin_info = [(bin_centers[i], bin_intensities[i]) for i in 1:n_bins]
        return feature_matrix, bin_info
    else
        @warn "No bins created after filtering"
        return nothing, nothing
    end
end

"""
    save_feature_matrix(feature_matrix::Matrix{Float64}, bin_info, output_dir::String) -> Tuple{String, String}

Saves the aggregated `2 x n_bins` feature matrix into two different CSV formats.

1.  **Simple Format (`feature_matrix_simple.csv`):** A two-column CSV with "mz" and "intensity".
2.  **Standard Format (`feature_matrix_standard.csv`):** A row-based format where m/z values are headers and there is a single data row for the aggregated spectrum.

# Arguments
- `feature_matrix::Matrix{Float64}`: The `2 x n_bins` matrix from `apply_peak_binning`.
- `bin_info`: The associated bin information (currently unused but kept for compatibility).
- `output_dir::String`: The directory where the output CSV files will be saved.

# Returns
- A tuple containing the paths to the two saved files.
"""
function save_feature_matrix(feature_matrix::Matrix{Float64}, bin_info, output_dir::String)
    # Save as simple CSV with m/z and intensity rows
    csv_path = joinpath(output_dir, "feature_matrix_simple.csv")
    
    open(csv_path, "w") do io
        write(io, "mz,intensity\n")
        for i in 1:size(feature_matrix, 2)
            mz = feature_matrix[1, i]
            intensity = feature_matrix[2, i]
            write(io, "$mz,$intensity\n")
        end
    end
    @info "Saved simple feature matrix: $csv_path"
    
    # Also save in a more standard format for MSI
    csv_path_standard = joinpath(output_dir, "feature_matrix_standard.csv")
    
    open(csv_path_standard, "w") do io
        write(io, "sample_type,")
        mz_headers = [@sprintf("mz_%.4f", feature_matrix[1, i]) for i in 1:size(feature_matrix, 2)]
        write(io, join(mz_headers, ",") * "\n")
        
        write(io, "aggregated_spectrum,")
        intensity_values = [feature_matrix[2, i] for i in 1:size(feature_matrix, 2)]
        write(io, join(string.(intensity_values), ",") * "\n")
    end
    @info "Saved standard format matrix: $csv_path_standard"
    
    return csv_path, csv_path_standard
end

function execute_full_preprocessing(spectra::Vector{MutableSpectrum}, params::Dict, 
                                   pipeline_steps::Vector{String}, reference_peaks::Dict, 
                                   mask_path::Union{String, Nothing}=nothing;
                                   progress_callback::Function=(step -> nothing))
    
    println("Starting preprocessing pipeline with $(length(spectra)) spectra")
    println("Steps: $(join(pipeline_steps, " -> "))")
    
    # These variables will be populated by the pipeline steps
    feature_matrix = nothing
    bin_definitions = nothing

    # Apply pipeline steps, modifying `spectra` in-place
    for step in pipeline_steps
        progress_callback(step)
        println("\n" * "-"^60)
        println("PROCESSING STEP: $step")
        println("-"^60)
        
        if step == "stabilization"
            println("  Applying intensity transformation (stabilization)")
            apply_intensity_transformation(spectra, get(params, :Stabilization, Dict()))
            
        elseif step == "baseline_correction"
            println("  Applying baseline correction")
            apply_baseline_correction(spectra, get(params, :BaselineCorrection, Dict()))
            
        elseif step == "smoothing"
            println("  Applying smoothing")
            apply_smoothing(spectra, get(params, :Smoothing, Dict()))
            
        elseif step == "peak_picking"
            println("  Applying peak picking")
            apply_peak_picking(spectra, get(params, :PeakPicking, Dict()))

        elseif step == "peak_selection"
            println("  Applying peak selection")
            apply_peak_selection(spectra, get(params, :PeakSelection, Dict()))
            
        elseif step == "calibration"
            println("  Applying calibration")
            apply_calibration(spectra, get(params, :Calibration, Dict()), reference_peaks)
            
        elseif step == "peak_alignment"
            println("  Applying peak alignment")
            apply_peak_alignment(spectra, get(params, :PeakAlignment, Dict()))
            
        elseif step == "normalization"
            println("  Applying normalization")
            apply_normalization(spectra, get(params, :Normalization, Dict()))
            
        elseif step == "peak_binning"
            println("  Applying peak binning")
            feature_matrix, bin_definitions = apply_peak_binning(spectra, get(params, :PeakBinning, Dict()))
            
        else
            @warn "Unknown step: $step, skipping"
        end
        
        println("✓ Completed step: $step")
    end
    
    return feature_matrix, bin_definitions
end