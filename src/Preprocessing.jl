# src/Preprocessing.jl
"""
This module provides a comprehensive workflow for mass spectrometry imaging (MSI) data
preprocessing, inspired by the functionality of the R package MALDIquant and Cardinal. 
It includes functions for quality control, intensity transformation, smoothing, 
baseline correction, normalization, peak picking, alignment, and feature matrix 
generation.
"""

# =============================================================================
# Dependencies
# =============================================================================

using Statistics  # For mean, median
using StatsBase   # For mad (Median Absolute Deviation)
using SavitzkyGolay # For SavitzkyGolay filtering
using Dates         # For now()
using CSV           # For writing CSV files
using DataFrames    # For creating dataframes
using ContinuousWavelets # For CWT peak detection
using ImageFiltering # For localmaxima in detect_peaks_wavelet
using Interpolations # For calibration
using Loess # For robust peak alignment

# =============================================================================
# Data Structures
# =============================================================================

"""
    FeatureMatrix

A struct to hold the final feature matrix generated from the preprocessing pipeline.
"""
struct FeatureMatrix
    matrix::Array{Float64,2}
    mz_bins::Vector{Tuple{Float64,Float64}}
    sample_ids::Vector{Int}
end

"""
A mutable struct to hold spectrum data. Using a mutable struct allows
in-place modification of fields (like intensity or m/z), which dramatically
reduces memory allocations compared to creating new immutable tuples at each step.
"""
mutable struct MutableSpectrum
    id::Int
    mz::AbstractVector{Float64}
    intensity::AbstractVector{Float64}
    peaks::Vector{NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), NTuple{6, Float64}}}
end

# --- Core Pipeline Structs ---

"""
An abstract type for all preprocessing steps. Allows for a modular and extensible
pipeline where users can define a sequence of operations.
"""
abstract type AbstractPreprocessingStep end

# --- Preprocessing Step Definitions ---

"""
    Calibration(; method=:internal_standards, ...)

A preprocessing step for mass calibration. Set a parameter to `nothing` to use an
auto-determined value from the data where applicable.
"""
struct Calibration <: AbstractPreprocessingStep
    method::Symbol
    internal_standards::Union{Dict{Float64, String}, Nothing}
    base_peak_mz_references::Union{Vector{Float64}, Nothing}
    ppm_tolerance::Union{Float64, Nothing}
    fit_order::Int

    function Calibration(; method=:internal_standards, internal_standards=nothing, base_peak_mz_references=nothing, ppm_tolerance=nothing, fit_order=2)
        new(method, internal_standards, base_peak_mz_references, ppm_tolerance, fit_order)
    end
end

"""
    BaselineCorrection(; method=:snip, ...)

A preprocessing step for baseline correction.
"""
struct BaselineCorrection <: AbstractPreprocessingStep
    method::Symbol
    iterations::Union{Int, Nothing} # For SNIP
    window::Union{Int, Nothing}     # For median

    function BaselineCorrection(; method=:snip, iterations=nothing, window=nothing)
        new(method, iterations, window)
    end
end

"""
    Smoothing(; method=:savitzky_golay, ...)

A preprocessing step for spectral smoothing.
"""
struct Smoothing <: AbstractPreprocessingStep
    method::Symbol
    window::Union{Int, Nothing}
    order::Union{Int, Nothing} # For Savitzky-Golay

    function Smoothing(; method=:savitzky_golay, window=nothing, order=nothing)
        new(method, window, order)
    end
end



"""
    Normalization(; method=:tic)

A preprocessing step for intensity normalization.
"""
struct Normalization <: AbstractPreprocessingStep
    method::Symbol

    function Normalization(; method=:tic)
        new(method)
    end
end

"""
    PeakPicking(; method=nothing, ...)

A preprocessing step for peak detection.
"""
struct PeakPicking <: AbstractPreprocessingStep
    method::Union{Symbol, Nothing} # :profile, :wavelet, :centroid
    snr_threshold::Union{Float64, Nothing}
    half_window::Union{Int, Nothing}
    min_peak_prominence::Union{Float64, Nothing}
    merge_peaks_tolerance::Union{Float64, Nothing}
    min_peak_width_ppm::Union{Float64, Nothing}
    max_peak_width_ppm::Union{Float64, Nothing}
    min_peak_shape_r2::Union{Float64, Nothing}

    function PeakPicking(; method=nothing, snr_threshold=nothing, half_window=nothing, min_peak_prominence=nothing, merge_peaks_tolerance=nothing, min_peak_width_ppm=nothing, max_peak_width_ppm=nothing, min_peak_shape_r2=nothing)
        new(method, snr_threshold, half_window, min_peak_prominence, merge_peaks_tolerance, min_peak_width_ppm, max_peak_width_ppm, min_peak_shape_r2)
    end
end

"""
    PeakAlignment(; method=:lowess, ...)

A preprocessing step for peak alignment.
"""
struct PeakAlignment <: AbstractPreprocessingStep
    method::Symbol
    span::Union{Float64, Nothing}
    tolerance::Union{Float64, Nothing}
    tolerance_unit::Union{Symbol, Nothing}
    max_shift_ppm::Union{Float64, Nothing}
    min_matched_peaks::Union{Int, Nothing}

    function PeakAlignment(; method=:lowess, span=nothing, tolerance=nothing, tolerance_unit=nothing, max_shift_ppm=nothing, min_matched_peaks=nothing)
        new(method, span, tolerance, tolerance_unit, max_shift_ppm, min_matched_peaks)
    end
end

"""
    PeakSelection(; frequency_threshold=nothing, ...)

A preprocessing step for peak selection (filtering).
"""
struct PeakSelection <: AbstractPreprocessingStep
    frequency_threshold::Union{Float64, Nothing}
    min_snr::Union{Float64, Nothing}
    min_fwhm_ppm::Union{Float64, Nothing}
    max_fwhm_ppm::Union{Float64, Nothing}
    min_shape_r2::Union{Float64, Nothing}
    correlation_threshold::Union{Float64, Nothing}

    function PeakSelection(; frequency_threshold=nothing, min_snr=nothing, min_fwhm_ppm=nothing, max_fwhm_ppm=nothing, min_shape_r2=nothing, correlation_threshold=nothing)
        new(frequency_threshold, min_snr, min_fwhm_ppm, max_fwhm_ppm, min_shape_r2, correlation_threshold)
    end
end

"""
    PeakBinningParams(; method=:adaptive, ...)

A preprocessing step for peak binning.
"""
struct PeakBinningParams <: AbstractPreprocessingStep
    method::Symbol
    tolerance::Union{Float64, Nothing}
    tolerance_unit::Union{Symbol, Nothing}
    frequency_threshold::Union{Float64, Nothing}
    min_peak_per_bin::Union{Int, Nothing}
    max_bin_width_ppm::Union{Float64, Nothing}
    intensity_weighted_centers::Bool
    num_uniform_bins::Union{Int, Nothing}

    function PeakBinningParams(; method=:adaptive, tolerance=nothing, tolerance_unit=nothing, frequency_threshold=nothing, min_peak_per_bin=nothing, max_bin_width_ppm=nothing, intensity_weighted_centers=true, num_uniform_bins=nothing)
        new(method, tolerance, tolerance_unit, frequency_threshold, min_peak_per_bin, max_bin_width_ppm, intensity_weighted_centers, num_uniform_bins)
    end
end

# =============================================================================
# 0) Quality Control (QC)
# =============================================================================

"""
    qc_is_empty(mz, intensity) -> Bool

Returns `true` if the spectrum is empty or contains no finite intensity data.
"""
qc_is_empty(mz::AbstractVector, intensity::AbstractVector)::Bool =
    isempty(mz) || isempty(intensity) || all(!isfinite, intensity)

"""
    qc_is_regular(mz) -> Bool

Checks that the m/z axis is monotonically non-decreasing.
"""
function qc_is_regular(mz::AbstractVector)
    n = length(mz)
    n < 2 && return true
    @inbounds for i in 2:n
        if mz[i] < mz[i-1]
            return false
        end
    end
    return true
end

"""
    validate_spectrum(mz, intensity) -> Bool

Performs comprehensive validation on a mass spectrum (m/z and intensity arrays).
Returns `true` if the spectrum is valid, `false` otherwise.

Checks include:
- Both `mz` and `intensity` are non-empty.
- `mz` and `intensity` have the same length.
- All `m/z` values are finite and non-negative.
- All `intensity` values are finite and non-negative.
- `m/z` values are strictly increasing (no duplicates or decreasing values).
"""
function validate_spectrum(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real})::Bool
    # 1. Check for empty vectors
    if isempty(mz) || isempty(intensity)
        @warn "Spectrum validation failed: m/z or intensity vector is empty."
        return false
    end

    # 2. Check for length mismatch
    if length(mz) != length(intensity)
        @warn "Spectrum validation failed: m/z and intensity vectors have different lengths ($(length(mz)) vs $(length(intensity)))"
        return false
    end

    # 3. Check for NaN/Inf in mz and non-negativity
    if any(!isfinite, mz) || any(<(0), mz)
        @warn "Spectrum validation failed: m/z vector contains NaN, Inf, or negative values."
        return false
    end

    # 4. Check for NaN/Inf in intensity and non-negativity
    if any(!isfinite, intensity) || any(<(0), intensity)
        #@warn "Spectrum validation failed: intensity vector contains NaN, Inf, or negative values."
        return false
    end

    # 5. Check for strictly increasing m/z values (and thus no duplicates)
    if !qc_is_regular(mz)
        #@warn "Spectrum validation failed: m/z vector is not strictly increasing (contains duplicates or decreasing values)."
        return false
    end

    return true
end

# =============================================================================
# 1) Intensity Transformation & Smoothing
# =============================================================================

"""
    transform_intensity(intensity; method=:sqrt) -> Vector

Applies a variance-stabilizing transformation to the intensity vector.
"""
function transform_intensity(intensity::AbstractVector{<:Real}; method::Symbol=:sqrt)
    if method === :sqrt
        return sqrt.(max.(zero(eltype(intensity)), intensity))
    elseif method === :log1p
        return log1p.(max.(zero(eltype(intensity)), intensity))
    elseif method === :log
        return log.(max.(eps(eltype(intensity)), intensity))
    elseif method === :log2
        return log2.(max.(eps(eltype(intensity)), intensity))
    elseif method === :log10
        return log10.(max.(eps(eltype(intensity)), intensity))
    else
        return collect(float.(intensity))
    end
end

"""
    smooth_spectrum(y::AbstractVector{<:Real}; method::Symbol=:savitzky_golay, window::Int=9, order::Int=2) -> Vector

Applies a smoothing filter to the intensity data.

# Arguments
- `y`: The intensity data.
- `method`: The smoothing method (:savitzky_golay or :moving_average).
- `window`: The window size for the filter.
- `order`: The polynomial order for Savitzky-Golay
"""
function smooth_spectrum(y::AbstractVector{<:Real}; method::Symbol=:savitzky_golay, window::Int=9, order::Int=2)
    if window < 3
        throw(ArgumentError("Window size must be at least 3"))
    end
    if order >= window
        throw(ArgumentError("Polynomial order must be less than window size"))
    end
    if order < 0
        throw(ArgumentError("Polynomial order cannot be negative"))
    end
    if method === :savitzky_golay
        win = isodd(window) ? window : window + 1
        if length(y) < win
            return y # Cannot smooth if data is smaller than window
        end
        res = SavitzkyGolay.savitzky_golay(collect(float.(y)), win, order)
        return res.y
    elseif method === :moving_average
        return moving_average_smooth(y, window)
    else
        @warn "Unsupported smoothing method: $method. Returning original intensity."
        return collect(float.(y))
    end
end

"""
    moving_average_smooth(y::AbstractVector{<:Real}, window::Int) -> Vector

Applies a simple moving average filter to the intensity data.
"""
function moving_average_smooth(y::AbstractVector{<:Real}, window::Int)
    n = length(y)
    if n < window
        return collect(float.(y)) # Cannot smooth if data is smaller than window
    end
    
    smoothed_y = zeros(Float64, n)
    half_window = div(window, 2)
    
    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        smoothed_y[i] = mean(@view y[start_idx:end_idx])
    end
    return smoothed_y
end

# =============================================================================
# 2) Baseline Correction
# =============================================================================

"""
    _snip_baseline_impl(y, iterations=100) -> Vector

Estimates the baseline of a spectrum using the SNIP algorithm (internal implementation).
"""
function _snip_baseline_impl(y::AbstractVector{<:Real}; iterations::Int=100)
    n = length(y)
    
    # Initialize two buffers. b1 holds the current baseline estimate, b2 for the next.
    # Always convert to Float64 to ensure type stability and avoid copying if already correct type
    b1 = collect(float.(y))
    b2 = similar(b1)        
    
    current_b = b1
    next_b = b2

    for k in 1:iterations
        # Calculate next baseline estimate into `next_b` based on `current_b`
        # Boundary conditions
        if n > 1
            next_b[1] = min(current_b[1], current_b[2])
            next_b[n] = min(current_b[n], current_b[n-1])
        end

        @inbounds for i in 2:n-1
            next_b[i] = min(current_b[i], 0.5 * (current_b[i-1] + current_b[i+1]))
        end
        
        # Swap references for the next iteration (no data copy here)
        current_b, next_b = next_b, current_b 
    end
    
    # Return the final baseline estimate (which is in current_b after the last swap)
    return current_b
end

"""
    convex_hull_baseline(y) -> Vector

Estimates the baseline of a spectrum using the convex hull algorithm.
"""
function convex_hull_baseline(y::AbstractVector{<:Real})
    n = length(y)
    if n < 3
        return zeros(Float64, n)
    end

    # Find upper convex hull points
    hull_indices = Int[1]
    for i in 2:n
        while length(hull_indices) >= 2 &&
              (y[hull_indices[end]] - y[hull_indices[end-1]]) * (i - hull_indices[end]) <= 
              (y[i] - y[hull_indices[end]]) * (hull_indices[end] - hull_indices[end-1])
            pop!(hull_indices)
        end
        push!(hull_indices, i)
    end

    # Interpolate between hull points
    baseline = zeros(Float64, n)
    for i in 1:(length(hull_indices)-1)
        idx1 = hull_indices[i]
        idx2 = hull_indices[i+1]
        
        # Linear interpolation
        for j in idx1:idx2
            baseline[j] = y[idx1] + (y[idx2] - y[idx1]) * (j - idx1) / (idx2 - idx1)
        end
    end
    return baseline
end

"""
    median_baseline(y; window=20) -> Vector

Estimates the baseline of a spectrum using a moving median filter.
"""
function median_baseline(y::AbstractVector{<:Real}; window::Int=20)
    n = length(y)
    if n < window
        return zeros(Float64, n)
    end

    baseline = zeros(Float64, n)
    half_window = div(window, 2)

    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        baseline[i] = median(@view y[start_idx:end_idx])
    end
    return baseline
end

"""
    apply_baseline_correction(y::AbstractVector{<:Real}; method::Symbol=:snip, iterations::Int=100, window::Int=20) -> Vector

Applies a baseline correction algorithm to the intensity data.

# Arguments
- `y`: The intensity data.
- `method`: The baseline correction method (:snip, :convex_hull, or :median).
- `iterations`: Iterations for SNIP method.
- `window`: Window size for median method.
"""
function apply_baseline_correction(y::AbstractVector{<:Real}; method::Symbol=:snip, iterations::Int=100, window::Int=20)
    if method === :snip
        return _snip_baseline_impl(y, iterations=iterations)
    elseif method === :convex_hull
        return convex_hull_baseline(y)
    elseif method === :median
        return median_baseline(y, window=window)
    else
        @warn "Unsupported baseline correction method: $method. Returning zero baseline."
        return zeros(Float64, length(y))
    end
end

# =============================================================================
# 3) Intensity Normalization
# =============================================================================

"""
    tic_normalize(y) -> Vector

Normalizes spectrum intensities to the Total Ion Current (TIC).
"""
function tic_normalize(y::AbstractVector{<:Real})
    s = sum(y)
    return s <= 0 ? collect(float.(y)) : collect(float.(y)) ./ s
end

"""
    pqn_normalize(M) -> Matrix

Performs Probabilistic Quotient Normalization (PQN) on a matrix of spectra.
"""
function pqn_normalize(M::AbstractMatrix{<:Real})
    M_float = collect(float.(M))
    ref = mapslices(median, M_float; dims=2)[:,1]
    Q = similar(M_float)
    @inbounds for j in axes(M_float, 2)
        Q[:, j] = M_float[:, j] ./ (ref .+ eps(eltype(M_float)))
    end
    s = [median(@view Q[:, j]) for j in axes(Q, 2)]
    @inbounds for j in axes(M_float, 2)
        M_float[:, j] ./= (s[j] + eps(eltype(M_float)))
    end
    return M_float
end

"""
    median_normalize(y) -> Vector

Normalizes spectrum intensities by dividing by the median intensity.
"""
function median_normalize(y::AbstractVector{<:Real})
    m = median(y)
    return m <= 0 ? collect(float.(y)) : collect(float.(y)) ./ m
end

"""
    rms_normalize(y) -> Vector

Normalizes spectrum intensities to the Root Mean Square (RMS).
"""
function rms_normalize(y::AbstractVector{<:Real})
    s = sqrt(sum(abs2, y) / length(y))
    return s <= 0 ? collect(float.(y)) : collect(float.(y)) ./ s
end

"""
    apply_normalization(y::AbstractVector{<:Real}; method::Symbol=:tic) -> Vector

Applies a per-spectrum normalization algorithm to the intensity data.

# Arguments
- `y`: The intensity data.
- `method`: The normalization method (:tic, :median, :rms, or :none).
"""
function apply_normalization(y::AbstractVector{<:Real}; method::Symbol=:tic)
    if method === :tic
        return tic_normalize(y)
    elseif method === :median
        return median_normalize(y)
    elseif method === :rms
        return rms_normalize(y)
    elseif method === :none
        return collect(float.(y))
    else
        @warn "Unsupported normalization method: $method. Returning original intensity."
        return collect(float.(y))
    end
end

"""
    remove_matrix_peaks_from_spectrum(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real},
                                      matrix_peak_mzs::AbstractVector{<:Real};
                                      tolerance::Float64=0.002, tolerance_unit::Symbol=:mz,
                                      removal_method::Symbol=:zero_out) -> Vector{Float64}

Removes or reduces intensity around specified matrix peaks in a single spectrum.

# Arguments
- `mz`: The m/z vector of the spectrum.
- `intensity`: The intensity vector of the spectrum.
- `matrix_peak_mzs`: A list of m/z values identified as matrix peaks.
- `tolerance`: The m/z tolerance for matching matrix peaks.
- `tolerance_unit`: Unit of tolerance (:mz or :ppm).
- `removal_method`: How to remove the peaks (:zero_out or :subtract).

# Returns
- `Vector{Float64}`: The modified intensity vector.
"""
function remove_matrix_peaks_from_spectrum(mz::AbstractVector{<:Real}, intensity::AbstractVector{<:Real},
                                           matrix_peak_mzs::AbstractVector{<:Real};
                                           tolerance::Float64=0.002, tolerance_unit::Symbol=:mz,
                                           removal_method::Symbol=:zero_out)
    modified_intensity = copy(intensity)
    
    for matrix_mz in matrix_peak_mzs
        # Calculate dynamic tolerance if in PPM
        current_tolerance = (tolerance_unit == :ppm) ? (matrix_mz * tolerance / 1e6) : tolerance
        
        # Find indices within the tolerance window
        indices_to_modify = findall(m -> abs(m - matrix_mz) <= current_tolerance, mz)
        
        if !isempty(indices_to_modify)
            if removal_method == :zero_out
                modified_intensity[indices_to_modify] .= 0.0
            elseif removal_method == :subtract
                # Baseline-aware subtraction: replace peak with a line connecting its "feet"
                idx_start = indices_to_modify[1]
                idx_end = indices_to_modify[end]
                
                # Ensure we are not at the very edge of the spectrum
                if idx_start > 1 && idx_end < length(modified_intensity)
                    y1 = modified_intensity[idx_start - 1]
                    y2 = modified_intensity[idx_end + 1]
                    x1 = idx_start - 1
                    x2 = idx_end + 1

                    # Linearly interpolate the baseline under the peak
                    for i in idx_start:idx_end
                        # y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
                        baseline_val = y1 + (y2 - y1) * (i - x1) / (x2 - x1)
                        # Set the intensity to the baseline, but don't increase it (e.g., if baseline is above signal)
                        modified_intensity[i] = min(modified_intensity[i], baseline_val)
                    end
                else
                    # If peak is at the edge, we can't interpolate, so just zero it out
                    modified_intensity[indices_to_modify] .= 0.0
                end
            else
                @warn "Unsupported matrix peak removal method: $removal_method. Skipping."
            end
        end
    end
    return modified_intensity
end

function identify_matrix_peaks_from_blanks(
    msi_data::MSIData,
    blank_spectrum_tag::Symbol,
    snr_threshold::Float64;
    frequency_threshold::Float64=0.5, # e.g., peak must be in 50% of blanks
    bin_tolerance::Float64=0.005, # m/z tolerance for binning blank peaks
    bin_tolerance_unit::Symbol=:mz
    )::Vector{Float64}
    if !(0 < frequency_threshold <= 1.0)
        throw(ArgumentError("`frequency_threshold` must be between 0 and 1 (exclusive of 0). Got: $frequency_threshold"))
    end
    println("Identifying matrix peaks from blank spectra (tag: $blank_spectrum_tag)...")
    blank_indices = Int[]
    for (i, meta) in enumerate(msi_data.spectra_metadata)
        if meta.type == blank_spectrum_tag
            push!(blank_indices, i)
        end
    end

    if isempty(blank_indices)
        @warn "No blank spectra found with tag: $blank_spectrum_tag. Cannot identify matrix peaks."
        return Float64[]
    end

    all_blank_peaks = Vector{NamedTuple}[]
    num_blank_spectra = 0

    # Collect peaks from each blank spectrum
    _iterate_spectra_fast(msi_data, blank_indices) do idx, mz, intensity
        num_blank_spectra += 1
        if !validate_spectrum(mz, intensity)
            @warn "Blank spectrum $idx is invalid, skipping peak detection for it."
            push!(all_blank_peaks, NamedTuple[]) # Add empty list to maintain count
            return
        end
        # Assuming profile mode for matrix peak detection
        peaks = detect_peaks_profile(mz, intensity; snr_threshold=snr_threshold)
        push!(all_blank_peaks, peaks)
    end

    if isempty(all_blank_peaks) || all(isempty, all_blank_peaks)
        @warn "No peaks detected in any blank spectra. Cannot identify matrix peaks."
        return Float64[]
    end

    # Flatten all peaks from blank spectra for initial binning
    flat_blank_peaks = NamedTuple[]
    for peaks_in_spec in all_blank_peaks
        append!(flat_blank_peaks, peaks_in_spec)
    end
    sort!(flat_blank_peaks, by=p->p.mz)

    # Bin the detected peaks from all blanks to find common m/z features
    # This is a simplified binning for matrix peak identification
    binned_matrix_features = Dict{Float64, Int}() # mz_center => count of spectra it appeared in

    i = 1
    while i <= length(flat_blank_peaks)
        current_bin_start_idx = i
        current_peak = flat_blank_peaks[i]
        
        # Calculate dynamic tolerance if in PPM
        current_bin_mz = current_peak.mz
        tol_val = (bin_tolerance_unit == :ppm) ? (current_bin_mz * bin_tolerance / 1e6) : bin_tolerance
        
        j = i + 1
        while j <= length(flat_blank_peaks) && (flat_blank_peaks[j].mz - current_peak.mz) <= tol_val
            j += 1
        end
        current_bin_end_idx = j - 1

        # Calculate a representative m/z for the bin (e.g., intensity-weighted average)
        peaks_in_bin = flat_blank_peaks[current_bin_start_idx:current_bin_end_idx]
        if !isempty(peaks_in_bin)
            sum_intensity = sum(p.intensity for p in peaks_in_bin)
            if sum_intensity > 0
                bin_center_mz = sum(p.mz * p.intensity for p in peaks_in_bin) / sum(sum_intensity)
            else
                bin_center_mz = mean(p.mz for p in peaks_in_bin)
            end
            
            # Check how many blank spectra this feature appeared in
            spectra_count = 0
            # For each blank spectrum, check if it contains a peak within the current bin's tolerance
            for spec_peaks in all_blank_peaks
                if any(p -> abs(p.mz - bin_center_mz) <= tol_val, spec_peaks)
                    spectra_count += 1
                end
            end
            binned_matrix_features[bin_center_mz] = spectra_count
        end
        i = j
    end

    # Filter based on frequency_threshold
    matrix_peak_mzs = Float64[]
    min_spectra_count = ceil(Int, num_blank_spectra * frequency_threshold)

    for (mz_center, count) in binned_matrix_features
        if count >= min_spectra_count
            push!(matrix_peak_mzs, mz_center)
        end
    end
    sort!(matrix_peak_mzs)

    println("Identified $(length(matrix_peak_mzs)) potential matrix peaks from $(num_blank_spectra) blank spectra.")
    return matrix_peak_mzs
end


# =============================================================================
# 4) Peak Detection
# =============================================================================

"""
    detect_peaks_profile(mz, y; ...) -> Vector{NamedTuple}

Enhanced peak detection for profile-mode spectra with advanced filtering and quality metrics.

Returns a vector of NamedTuples, each representing a detected peak with:
- `mz`: m/z value of the peak
- `intensity`: Intensity of the peak
- `fwhm`: Full Width at Half Maximum (Δm)
- `shape_r2`: Pseudo R^2 for peak shape (simplified)
- `snr`: Signal-to-Noise Ratio
- `prominence`: Peak prominence
"""
function detect_peaks_profile(mz::AbstractVector{<:Real}, y::AbstractVector{<:Real}; 
                              half_window::Int=10,
                              snr_threshold::Float64=2.0,
                              min_peak_prominence::Float64=0.1,
                              merge_peaks_tolerance::Float64=0.002,
                              # New parameters for quality filtering (used for calculation, not filtering here)
                              min_peak_width_ppm::Float64=0.0, # Not used for filtering in this function
                              max_peak_width_ppm::Float64=Inf, # Not used for filtering in this function
                              min_peak_shape_r2::Float64=0.0 # Not used for filtering in this function
                              )
    n = length(y)
    n < 3 && return NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]
    
    noise_level = mad(y, normalize=true) + eps(Float64)
    ys = smooth_spectrum(y; method=:savitzky_golay, window=max(5, 2*half_window+1), order=2) # Use smoothed data for detection

    candidate_peak_indices = Int[]
    for i in 2:n-1
        left = max(1, i - half_window)
        right = min(n, i + half_window)
        
        # Prominence check
        prominence = ys[i] - max(minimum(@view ys[left:i]), minimum(@view ys[i:right]))
        
        if ys[i] >= maximum(@view ys[left:right]) && 
           (ys[i] > snr_threshold * noise_level) &&
           (prominence > min_peak_prominence * ys[i])
            push!(candidate_peak_indices, i)
        end
    end
    
    # Merge close peaks
    if !isempty(candidate_peak_indices) && merge_peaks_tolerance > 0
        merged_indices = [candidate_peak_indices[1]]
        for i in 2:length(candidate_peak_indices)
            if (mz[candidate_peak_indices[i]] - mz[last(merged_indices)]) > merge_peaks_tolerance
                push!(merged_indices, candidate_peak_indices[i])
            elseif y[candidate_peak_indices[i]] > y[last(merged_indices)]
                merged_indices[end] = candidate_peak_indices[i] # Replace with more intense peak
            end
        end
        candidate_peak_indices = merged_indices
    end

    detected_peaks = NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]
    for p_idx in candidate_peak_indices
        peak_mz = float(mz[p_idx])
        peak_int = float(y[p_idx])
        
        fwhm_delta_m = calculate_robust_fwhm(mz, y, p_idx)
        fwhm_ppm = if isnan(fwhm_delta_m) || fwhm_delta_m <= 0
            0.0
        else
            1e6 * fwhm_delta_m / peak_mz
        end
        shape_r2 = _fit_gaussian_and_r2(mz, y, p_idx, half_window)
        
        peak_snr = peak_int / noise_level # Recalculate SNR based on final peak_int
        
        left = max(1, p_idx - half_window)
        right = min(n, p_idx + half_window)
        prominence = ys[p_idx] - max(minimum(@view ys[left:p_idx]), minimum(@view ys[p_idx:right]))

        push!(detected_peaks, (mz=peak_mz, intensity=peak_int, fwhm=fwhm_ppm, shape_r2=shape_r2, snr=peak_snr, prominence=prominence))
    end
    
    return detected_peaks
end

"""
    detect_peaks_wavelet(mz, intensity; ...) -> Vector{NamedTuple}

Detects peaks using Continuous Wavelet Transform (CWT).

Returns a vector of NamedTuples, each representing a detected peak with:
- `mz`: m/z value of the peak
- `intensity`: Intensity of the peak
- `snr`: Signal-to-Noise Ratio (simplified)
"""
function detect_peaks_wavelet(mz::AbstractVector, intensity::AbstractVector; scales=1:10, snr_threshold=3.0, half_window=10)::Vector{NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}}
    n = length(intensity)
    n < 10 && return NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]
    
    # Compute CWT
    cwt_res = ContinuousWavelets.cwt(intensity, ContinuousWavelets.morl) # Morlet is good for peaks
    
    noise_level = mad(intensity, normalize=true) + eps(Float64)
    
    candidate_indices = Set{Int}() # Use a Set to store unique peak indices

    # Find local maxima in the CWT coefficients (magnitude)
    # This assumes `localmaxima` is available (e.g., from ImageFiltering, often a dependency of ContinuousWavelets).
    abs_cwt_res = abs.(cwt_res)
    # Iterate over all (m/z index, scale index) pairs that are local maxima in the CWT matrix
    for (m_idx, scale_idx) in Tuple.(localmaxima(abs_cwt_res))
        # Ensure m_idx is not at the very edges of the intensity array to avoid index out of bounds
        if m_idx > 1 && m_idx < n
            # Check if this CWT maximum corresponds to a local maximum in the original intensity
            # AND if the original intensity is above the SNR threshold
            if intensity[m_idx] > intensity[m_idx-1] && intensity[m_idx] > intensity[m_idx+1] && intensity[m_idx] > snr_threshold * noise_level
                push!(candidate_indices, m_idx)
            end
        end
    end
    peak_indices = collect(candidate_indices)
    sort!(peak_indices) # Ensure order

    detected_peaks = NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]
    for p_idx in peak_indices
        peak_mz = float(mz[p_idx])
        peak_int = float(intensity[p_idx])
        peak_snr = peak_int / noise_level
        
        fwhm_delta_m = calculate_robust_fwhm(mz, intensity, p_idx)
        fwhm_ppm = if isnan(fwhm_delta_m) || fwhm_delta_m <= 0
            0.0
        else
            1e6 * fwhm_delta_m / peak_mz
        end
        shape_r2 = _fit_gaussian_and_r2(mz, intensity, p_idx, half_window)

        push!(detected_peaks, (mz=peak_mz, intensity=peak_int, fwhm=fwhm_ppm, shape_r2=shape_r2, snr=peak_snr, prominence=peak_int))
    end
    
    return detected_peaks
end

"""
    detect_peaks_centroid(mz, y; ...) -> Vector{NamedTuple}

Filters peaks in centroid-mode data based on intensity threshold.

Returns a vector of NamedTuples, each representing a detected peak with:
- `mz`: m/z value of the peak
- `intensity`: Intensity of the peak
"""
function detect_peaks_centroid(mz::AbstractVector{<:Real}, y::AbstractVector{<:Real}; snr_threshold::Float64=0.0)
    noise_level = mad(y, normalize=true) + eps(Float64)
    detected_peaks = NamedTuple{(:mz, :intensity, :fwhm, :shape_r2, :snr, :prominence), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]
    
    for i in eachindex(mz)
        snr = noise_level > 0 ? y[i] / noise_level : y[i] > 0 ? Inf : 0.0
        if snr >= snr_threshold
            # For centroid data, FWHM and shape are not applicable. Return placeholders.
            push!(detected_peaks, (mz=float(mz[i]), intensity=float(y[i]), fwhm=0.0, shape_r2=1.0, snr=snr, prominence=y[i]))
        end
    end
    return detected_peaks
end

# =============================================================================
# 5) Peak Alignment & Calibration
# =============================================================================

"""
    align_peaks_lowess(ref_mz, tgt_mz; ...)

Enhanced peak alignment with PPM tolerance and other constraints.
"""
function align_peaks_lowess(ref_mz::Vector{<:Real}, tgt_mz::Vector{<:Real};
                            method::Symbol=:linear, # :linear, :lowess, or :ransac
                            span::Float64=0.75, # Span for LOWESS
                            tolerance::Float64=0.002,
                            tolerance_unit::Symbol=:mz,
                            max_shift_ppm::Float64=50.0,
                            min_matched_peaks::Int=5)
    pairs = Tuple{Float64,Float64}[]
    i = 1; j = 1
    while i <= length(tgt_mz) && j <= length(ref_mz)
        # Dynamic tolerance for PPM
        tol = (tolerance_unit == :ppm) ? (ref_mz[j] * tolerance / 1e6) : tolerance
        dt = tgt_mz[i] - ref_mz[j]
        
        if abs(dt) <= tol
            # Max shift check to avoid spurious matches
            if abs(dt) * 1e6 / ref_mz[j] <= max_shift_ppm
                push!(pairs, (float(tgt_mz[i]), float(ref_mz[j])))
            end
            i += 1; j += 1
        elseif dt < 0
            i += 1
        else
            j += 1
        end
    end

    if length(pairs) < min_matched_peaks
        #@warn "Too few matching peaks ($(length(pairs)) < $min_matched_peaks). Returning identity function."
        return x -> float.(x)
    end

    t = [p[1] for p in pairs] # Target m/z (x-axis)
    r = [p[2] for p in pairs] # Reference m/z (y-axis) 
    
    if method == :linear
        itp = linear_interpolation(t, r, extrapolation_bc=Line())
        return itp
    elseif method == :lowess
        try
            model = loess(t, r; span=span)
            # Predict over the original target m/z values (t)
            predicted_r = Loess.predict(model, t)
            # Create a linear interpolation from original t and Loess-predicted r,
            # allowing linear extrapolation.
            itp = linear_interpolation(t, predicted_r, extrapolation_bc=Line())
            return itp
        catch e
            @warn "LOWESS fitting failed: $e. Falling back to linear interpolation."
            itp = linear_interpolation(t, r, extrapolation_bc=Line())
            return itp
        end
    elseif method == :ransac
        # RANSAC implementation for linear model
        num_iterations = 100
        best_model_inliers_count = -1
        best_inliers_t = Float64[]
        best_inliers_r = Float64[]

        n_pairs = length(pairs)
        if n_pairs < 2 # Need at least 2 points to fit a line
            @warn "Not enough matched peaks ($n_pairs) for RANSAC. Falling back to linear interpolation."
            itp = linear_interpolation(t, r, extrapolation_bc=Line())
            return itp
        end

        for iter in 1:num_iterations
            # 1. Randomly select 2 points
            sample_indices = StatsBase.sample(1:n_pairs, 2, replace=false)
            p1_x, p1_y = pairs[sample_indices[1]]
            p2_x, p2_y = pairs[sample_indices[2]]

            # Avoid vertical line for simplicity in this basic implementation
            if isapprox(p1_x, p2_x, atol=1e-9)
                continue 
            end

            # 2. Fit a line y = mx + b
            m = (p2_y - p1_y) / (p2_x - p1_x)
            b = p1_y - m * p1_x

            current_inliers_t = Float64[]
            current_inliers_r = Float64[]
            current_inliers_count = 0

            # 3. Find inliers
            for (px, py) in pairs
                predicted_y = m * px + b
                residual = abs(py - predicted_y)
                
                # Determine inlier threshold based on tolerance_unit
                current_inlier_threshold = (tolerance_unit == :ppm) ? (px * tolerance / 1e6) : tolerance

                if residual <= current_inlier_threshold
                    current_inliers_count += 1
                    push!(current_inliers_t, px)
                    push!(current_inliers_r, py)
                end
            end

            # 4. Evaluate model
            if current_inliers_count > best_model_inliers_count
                best_model_inliers_count = current_inliers_count
                best_inliers_t = current_inliers_t
                best_inliers_r = current_inliers_r
            end
        end

        if best_model_inliers_count >= 2 # At least 2 inliers to form a line
            return linear_interpolation(best_inliers_t, best_inliers_r, extrapolation_bc=Line())
        else
            @warn "RANSAC failed to find a robust model (found $(best_model_inliers_count) inliers). Falling back to linear interpolation."
            itp = linear_interpolation(t, r, extrapolation_bc=Line())
            return itp
        end
    else
        @warn "Unsupported alignment method: $method. Falling back to linear interpolation."
        itp = linear_interpolation(t, r, extrapolation_bc=Line())
        return itp
    end
end

"""
    find_calibration_peaks(mz, intensity, reference_masses; ...)

Finds peaks that match a list of reference masses.
"""
function find_calibration_peaks(mz::AbstractVector, intensity::AbstractVector, reference_masses::AbstractVector; ppm_tolerance=20.0)
    matched_peaks = Dict{Float64, Float64}()
    
    # detect_peaks_profile returns Vector{NamedTuple}, so we need to extract mz values
    detected_peaks_list = detect_peaks_profile(mz, intensity)
    
    # Extract only the m/z values into a new vector for easier processing
    detected_mz_values = [p.mz for p in detected_peaks_list]
    
    for ref_mass in reference_masses
        tol = ref_mass * ppm_tolerance / 1e6
        # Find candidates within the detected m/z values
        candidates = findall(m -> abs(m - ref_mass) <= tol, detected_mz_values)
        if !isempty(candidates)
            # Find the closest detected peak to the reference mass among candidates
            closest_peak_idx = argmin(abs.(detected_mz_values[candidates] .- ref_mass))
            matched_peaks[ref_mass] = detected_mz_values[candidates[closest_peak_idx]]
        end
    end
    return matched_peaks
end

"""
    calibrate_spectra(spectra, internal_standards; ...)

Calibrates spectra using internal standards.
"""
function calibrate_spectra(spectra::Vector, internal_standards::Vector; ppm_tolerance=20.0)
    calibrated_spectra = similar(spectra)
    for (i, spec) in enumerate(spectra)
        mz, intensity = spec[1], spec[2]

        matched_peaks = find_calibration_peaks(mz, intensity, internal_standards; ppm_tolerance=ppm_tolerance)
        if length(matched_peaks) < 2
            @warn "Spectrum $i: Not enough calibration peaks found. Skipping."
            calibrated_spectra[i] = spec
            continue
        end
        measured = sort(collect(values(matched_peaks)))
        theoretical = sort(collect(keys(matched_peaks)))
        itp = linear_interpolation(measured, theoretical, extrapolation_bc=Line())
        
        new_mz = itp(mz)

        if length(spec) == 3
            calibrated_spectra[i] = (new_mz, intensity, spec[3])
        else
            calibrated_spectra[i] = (new_mz, intensity)
        end
    end
    return calibrated_spectra
end

# =============================================================================
# 6) Peak Binning & Feature Matrix Generation
# =============================================================================

"""
    bin_peaks(all_pk_mz::Vector{Vector{Float64}},
              all_pk_int::Vector{Vector{Float64}},
              params::PeakBinningParams) -> Tuple{FeatureMatrix, Vector{Tuple{Float64,Float64}}}

Enhanced peak binning with adaptive and PPM-based parameters, or uniform binning.

# Arguments
- `all_pk_mz`: A vector of m/z vectors for all spectra.
- `all_pk_int`: A vector of intensity vectors for all spectra.
- `params`: A `PeakBinningParams` struct.

# Returns
- `Tuple{FeatureMatrix, Vector{Tuple{Float64,Float64}}}`: A tuple containing the generated FeatureMatrix and the bin definitions.

# Thread Safety
The use of `Threads.@threads` in the `:adaptive` and `:uniform` methods is safe.
- In the `:adaptive` method, the loop is over the bins (`j` index). Each thread writes only to its assigned column `X[:, j]`, so there are no write conflicts between threads.
- In the `:uniform` method, the loop is over the spectra (`s_idx`). Writes to `X[s_idx, bin_idx]` could theoretically conflict if different peaks from the same spectrum (`s_idx`) are processed by different threads. However, the loop is over `s_idx`, meaning each thread handles a distinct spectrum, making writes to `X[s_idx, :]` exclusive to that thread and thus safe.
"""
function bin_peaks(spectra::Vector{MutableSpectrum},
                   params::PeakBinningParams)
    ns = length(spectra) # Number of spectra
    ns == 0 && return FeatureMatrix(zeros(0,0), Tuple{Float64,Float64}[], Int[]), Tuple{Float64,Float64}[]

    if params.method == :uniform
        println("Performing uniform binning with $(params.num_uniform_bins) bins.")
        # Determine global m/z range from all peaks
        min_mz = Inf
        max_mz = -Inf
        for s in spectra
            if !isempty(s.peaks)
                # Extract m/z values from peaks NamedTuple
                pk_mz_vec = [p.mz for p in s.peaks] 
                min_mz = min(min_mz, minimum(pk_mz_vec))
                max_mz = max(max_mz, maximum(pk_mz_vec))
            end
        end

        if !isfinite(min_mz) || !isfinite(max_mz) || min_mz == max_mz
            @warn "Could not determine valid m/z range for uniform binning. Returning empty FeatureMatrix."
            return FeatureMatrix(zeros(0,0), Tuple{Float64,Float64}[], Int[]), Tuple{Float64,Float64}[]
        end

        # Create uniform bins
        mz_edges = collect(range(min_mz, stop=max_mz, length=params.num_uniform_bins + 1))
        bin_definitions = Vector{Tuple{Float64,Float64}}(undef, params.num_uniform_bins)
        for i in 1:params.num_uniform_bins
            bin_definitions[i] = (mz_edges[i], mz_edges[i+1])
        end
        
        X = zeros(Float64, ns, params.num_uniform_bins)
        
        Threads.@threads for s_idx in 1:ns
            s = spectra[s_idx] # Get the current MutableSpectrum
            for p in s.peaks # Iterate over peaks directly
                m = p.mz
                i = p.intensity
                # Find which bin this peak belongs to
                bin_idx = searchsortedlast(mz_edges, m)
                if bin_idx > 0 && bin_idx <= params.num_uniform_bins
                    X[s_idx, bin_idx] = max(X[s_idx, bin_idx], i) # Take max intensity in bin
                end
            end
        end

        if params.frequency_threshold !== nothing && params.frequency_threshold > 0
            present_count = vec(sum(X .> 0, dims=1))
            min_count = ceil(Int, params.frequency_threshold * ns)
            keep_mask = findall(present_count .>= min_count)
            X = X[:, keep_mask]
            bin_definitions = bin_definitions[keep_mask]
        end

        return FeatureMatrix(X, bin_definitions, collect(1:ns)), bin_definitions

    elseif params.method == :adaptive
        println("Performing adaptive binning with tolerance $(params.tolerance) $(params.tolerance_unit).")
        # Collect all peaks with their intensities and original spectrum ID directly from MutableSpectrum objects
        all_peaks = Vector{Tuple{Float64, Float64, Int}}()
        sizehint!(all_peaks, sum(length(s.peaks) for s in spectra)) # Pre-allocate memory

        for s in spectra
            for p in s.peaks
                push!(all_peaks, (p.mz, p.intensity, s.id))
            end
        end
        sort!(all_peaks, by=p->p[1])
        
        isempty(all_peaks) && return FeatureMatrix(zeros(0,0), Tuple{Float64,Float64}[], Int[]), Tuple{Float64,Float64}[]

        # Create bins using indices into all_peaks to avoid copying large vectors
        bins_indices = Vector{UnitRange{Int}}() # Stores UnitRange (start_idx:end_idx) for each bin
        if isempty(all_peaks)
            # Handle empty all_peaks case if it was not caught earlier
            return FeatureMatrix(zeros(0,0), Tuple{Float64,Float64}[], Int[]), Tuple{Float64,Float64}[]
        end

        current_bin_start_idx = 1
        for i in 2:length(all_peaks)
            p_prev = all_peaks[i-1]
            p_current = all_peaks[i]

            # Approximate bin center for tolerance calculation. Using current_bin_start_idx or p_current.mz is sufficient.
            # A more precise bin_center would require iterating over current_bin_peaks, which we are trying to avoid for memory.
            # For tolerance calculation, a single mz value (e.g., p_current.mz or all_peaks[current_bin_start_idx].mz) is often sufficient
            # because the m/z range within a small bin is typically very narrow.
            approx_bin_mz = p_prev[1] # Use the m/z of the previous peak in the bin as a proxy for bin center
            tol = (params.tolerance_unit == :ppm) ? (approx_bin_mz * params.tolerance / 1e6) : params.tolerance
            
            if (p_current[1] - p_prev[1]) <= tol
                # Peak is within tolerance, continue current bin
            else
                # Peak is outside tolerance, close current bin and start new one
                current_bin_end_idx = i - 1
                # Apply min_peak_per_bin filter when closing the bin
                if params.min_peak_per_bin === nothing || (current_bin_end_idx - current_bin_start_idx + 1) >= params.min_peak_per_bin
                    push!(bins_indices, current_bin_start_idx:current_bin_end_idx)
                end
                current_bin_start_idx = i
            end
        end
        # Add the last bin
        current_bin_end_idx = length(all_peaks)
        if params.min_peak_per_bin === nothing || (current_bin_end_idx - current_bin_start_idx + 1) >= params.min_peak_per_bin
            push!(bins_indices, current_bin_start_idx:current_bin_end_idx)
        end

        # Filter bins by max width
        if params.max_bin_width_ppm !== nothing
            filter!(range_idx -> begin
                bin_peaks_view = @view all_peaks[range_idx]
                min_mz_bin = first(bin_peaks_view)[1]
                max_mz_bin = last(bin_peaks_view)[1]
                bin_center_approx = (min_mz_bin + max_mz_bin) / 2
                (max_mz_bin - min_mz_bin) * 1e6 / bin_center_approx <= params.max_bin_width_ppm
            end, bins_indices)
        end

        # Create feature matrix
        X = zeros(Float64, ns, length(bins_indices))
        final_bins_boundaries = Vector{Tuple{Float64,Float64}}(undef, length(bins_indices))
        
        Threads.@threads for j in 1:length(bins_indices)
            range_idx = bins_indices[j]
            bin_peaks_view = @view all_peaks[range_idx] # Create a view into all_peaks for the current bin

            # Calculate bin center
            local bin_center
            if params.intensity_weighted_centers && sum(p[2] for p in bin_peaks_view) > 0
                weights = [p[2] for p in bin_peaks_view]
                bin_center = sum(p[1] * p[2] for p in bin_peaks_view) / sum(weights)
            else
                bin_center = mean(p[1] for p in bin_peaks_view)
            end
            
            final_bins_boundaries[j] = (first(bin_peaks_view)[1], last(bin_peaks_view)[1])

            for p in bin_peaks_view
                s_idx = p[3]
                # This is thread-safe because each thread writes to a unique column j
                X[s_idx, j] = max(X[s_idx, j], p[2])
            end
        end

        # Filter by frequency
        if params.frequency_threshold !== nothing && params.frequency_threshold > 0
            present_count = vec(sum(X .> 0, dims=1))
            min_count = ceil(Int, params.frequency_threshold * ns)
            keep_mask = findall(present_count .>= min_count)
            X = X[:, keep_mask]
            final_bins_boundaries = final_bins_boundaries[keep_mask]
        end

        return FeatureMatrix(X, final_bins_boundaries, collect(1:ns)), final_bins_boundaries
    else
        @warn "Unsupported binning method: $(params.method). Returning empty FeatureMatrix."
        return FeatureMatrix(zeros(0,0), Tuple{Float64,Float64}[], Int[]), Tuple{Float64,Float64}[]
    end
end
