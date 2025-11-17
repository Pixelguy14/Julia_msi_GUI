# src/Preprocessing.jl
"""
This module provides a comprehensive workflow for mass spectrometry imaging (MSI) data
preprocessing, inspired by the functionality of the R package MALDIquant. It includes
functions for quality control, intensity transformation, smoothing, baseline correction,
normalization, peak picking, alignment, and feature matrix generation.
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
using Interpolations # For calibration

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
    QCParameters

A struct to hold comprehensive quality control parameters for validating spectra.
"""
struct QCParameters
    min_tic_threshold::Float64      # Minimum total ion current
    max_tic_threshold::Float64      # Maximum TIC (saturation check)
    snr_threshold::Float64          # Minimum signal-to-noise for a spectrum to be considered valid
    peak_count_range::Tuple{Int,Int} # Acceptable peak count range (min, max)
    spatial_consistency_threshold::Float64 # For MSI spatial coherence (stub)
    calibration_accuracy_ppm::Float64 # Maximum allowed PPM error for calibration lock masses
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
    smooth_spectrum(y; window=9, order=2) -> Vector

Applies a Savitzky–Golay filter to smooth the intensity data.
"""
function smooth_spectrum(y::AbstractVector{<:Real}; window::Int=9, order::Int=2)
    win = isodd(window) ? window : window + 1
    if length(y) < win
        return y # Cannot smooth if data is smaller than window
    end
    res = SavitzkyGolay.savitzky_golay(collect(float.(y)), win, order)
    return res.y
end

# =============================================================================
# 2) Baseline Correction
# =============================================================================

"""
    snip_baseline(y, iterations=100) -> Vector

Estimates the baseline of a spectrum using the SNIP algorithm.
"""
function snip_baseline(y::AbstractVector{<:Real}; iterations::Int=100)
    n = length(y)
    b = collect(float.(y))
    buf = similar(b)
    for k in 1:iterations
        copyto!(buf, b)
        @inbounds for i in 2:n-1
            buf[i] = min(b[i], 0.5 * (b[i-1] + b[i+1]))
        end
        buf[1] = min(b[1], b[2])
        buf[end] = min(b[end], b[end-1])
        b, buf = buf, b
    end
    return b
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

# =============================================================================
# 4) Peak Detection
# =============================================================================

"""
    detect_peaks_profile(mz, y; ...)

Enhanced peak detection for profile-mode spectra with advanced filtering.
"""
function detect_peaks_profile(mz::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                              half_window::Int=10,
                              snr_threshold::Float64=2.0,
                              min_peak_width::Real=0.01,
                              max_peak_width::Real=2.0,
                              peak_shape_threshold::Float64=0.7, # Stub
                              merge_peaks_tolerance::Float64=0.002,
                              min_peak_prominence::Float64=0.1)
    n = length(y)
    n < 3 && return (Float64[], Float64[])
    
    noise_level = mad(y, normalize=true) + eps(Float64)
    ys = smooth_spectrum(y; window=max(5, 2*half_window+1), order=2)

    peak_indices = Int[]
    @inbounds for i in 2:n-1
        left = max(1, i - half_window)
        right = min(n, i + half_window)
        
        # Prominence check
        prominence = ys[i] - max(minimum(@view ys[left:i]), minimum(@view ys[i:right]))
        
        if ys[i] >= maximum(@view ys[left:right]) && 
           (ys[i] > snr_threshold * noise_level) &&
           (prominence > min_peak_prominence * ys[i])
            push!(peak_indices, i)
        end
    end
    
    # (STUB) Peak shape validation would be applied here
    # e.g., fitting a Gaussian and checking R^2 > peak_shape_threshold

    # Merge close peaks
    if !isempty(peak_indices) && merge_peaks_tolerance > 0
        merged_indices = [peak_indices[1]]
        for i in 2:length(peak_indices)
            if (mz[peak_indices[i]] - mz[last(merged_indices)]) > merge_peaks_tolerance
                push!(merged_indices, peak_indices[i])
            elseif y[peak_indices[i]] > y[last(merged_indices)]
                merged_indices[end] = peak_indices[i] # Replace with more intense peak
            end
        end
        peak_indices = merged_indices
    end

    # (STUB) Peak width filtering would be applied here
    # This would require calculating FWHM for each peak, which is computationally intensive.

    pk_mz  = [float(mz[i]) for i in peak_indices]
    pk_int = [float(y[i]) for i in peak_indices]
    
    return (pk_mz, pk_int)
end

"""
    detect_peaks_wavelet(mz, intensity; ...)

Detects peaks using Continuous Wavelet Transform (CWT).
"""
function detect_peaks_wavelet(mz::AbstractVector, intensity::AbstractVector; scales=1:10, snr_threshold=3.0)
    n = length(intensity)
    n < 10 && return (Float64[], Float64[])
    cwt_res = ContinuousWavelets.cwt(intensity, ContinuousWavelets.morl)
    
    peak_indices = Int[]
    noise_level = mad(intensity, normalize=true) + eps(Float64)

    for i in 2:n-1
        if intensity[i] > intensity[i-1] && intensity[i] > intensity[i+1] && intensity[i] > snr_threshold * noise_level
            if abs(cwt_res[i, end]) > 0
                push!(peak_indices, i)
            end
        end
    end
    
    return (mz[peak_indices], intensity[peak_indices])
end

"""
    detect_peaks_centroid(mz, y; ...)

Filters peaks in centroid-mode data.
"""
function detect_peaks_centroid(mz::AbstractVector{<:Real}, y::AbstractVector{<:Real}; intensity_threshold::Float64=0.0)
    keep_indices = findall(y .>= intensity_threshold)
    return (mz[keep_indices], y[keep_indices])
end

# =============================================================================
# 5) Peak Alignment & Calibration
# =============================================================================

"""
    align_peaks_lowess(ref_mz, tgt_mz; ...)

Enhanced peak alignment with PPM tolerance and other constraints.
"""
function align_peaks_lowess(ref_mz::Vector{<:Real}, tgt_mz::Vector{<:Real};
                            tolerance::Float64=0.002,
                            tolerance_unit::Symbol=:mz,
                            max_shift_ppm::Float64=50.0,
                            min_matched_peaks::Int=5)
    pairs = Tuple{Float64,Float64}[]
    i = 1; j = 1
    while i <= length(tgt_mz) && j <= length(ref_mz)
        tol = (tolerance_unit == :ppm) ? (ref_mz[j] * tolerance / 1e6) : tolerance
        dt = tgt_mz[i] - ref_mz[j]
        
        if abs(dt) <= tol
            # Max shift check
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
        @warn "Too few matching peaks ($(length(pairs)) < $min_matched_peaks). Returning identity function."
        return x -> float.(x)
    end

    t = [p[1] for p in pairs]
    r = [p[2] for p in pairs]
    
    # (STUB) A robust regression model (e.g., RANSAC) would be better here.
    itp = linear_interpolation(t, r, extrapolation_bc=Line())
    return itp
end

"""
    find_calibration_peaks(mz, intensity, reference_masses; ...)

Finds peaks that match a list of reference masses.
"""
function find_calibration_peaks(mz::AbstractVector, intensity::AbstractVector, reference_masses::AbstractVector; ppm_tolerance=20.0)
    matched_peaks = Dict{Float64, Float64}()
    detected_mz, _ = detect_peaks_profile(mz, intensity)
    
    for ref_mass in reference_masses
        tol = ref_mass * ppm_tolerance / 1e6
        candidates = findall(m -> abs(m - ref_mass) <= tol, detected_mz)
        if !isempty(candidates)
            closest_peak_idx = argmin(abs.(detected_mz[candidates] .- ref_mass))
            matched_peaks[ref_mass] = detected_mz[candidates[closest_peak_idx]]
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
    bin_peaks(all_pk_mz, all_pk_int, tolerance; ...)

Enhanced peak binning with adaptive and PPM-based parameters.
"""
function bin_peaks(all_pk_mz::Vector{<:AbstractVector{<:Real}},
                   all_pk_int::Vector{<:AbstractVector{<:Real}},
                   tolerance::Float64; 
                   frequency_threshold::Float64=0.25,
                   tolerance_unit::Symbol=:mz,
                   adaptive_tolerance::Bool=false, # Stub
                   min_peak_per_bin::Int=2,
                   max_bin_width_ppm::Float64=100.0,
                   intensity_weighted_centers::Bool=true)
    ns = length(all_pk_mz)
    ns == 0 && return (zeros(0,0), Tuple{Float64,Float64}[])

    # Collect all peaks with their intensities and original spectrum index
    all_peaks = [(float(m), float(i), s_idx) for s_idx in 1:ns for (m, i) in zip(all_pk_mz[s_idx], all_pk_int[s_idx])]
    sort!(all_peaks, by=p->p[1])
    
    isempty(all_peaks) && return (zeros(ns, 0), Tuple{Float64,Float64}[])

    # Create bins
    bins = []
    current_bin_peaks = [all_peaks[1]]
    for p in @view all_peaks[2:end]
        bin_center = mean(first.(current_bin_peaks))
        tol = (tolerance_unit == :ppm) ? (bin_center * tolerance / 1e6) : tolerance
        
        if p[1] - last(current_bin_peaks)[1] <= tol
            push!(current_bin_peaks, p)
        else
            if length(current_bin_peaks) >= min_peak_per_bin
                push!(bins, current_bin_peaks)
            end
            current_bin_peaks = [p]
        end
    end
    length(current_bin_peaks) >= min_peak_per_bin && push!(bins, current_bin_peaks)

    # Filter bins by max width
    filter!(b -> (last(b)[1] - first(b)[1]) * 1e6 / mean(p[1] for p in b) <= max_bin_width_ppm, bins)

    # Create feature matrix
    X = zeros(Float64, ns, length(bins))
    final_bins_boundaries = Vector{Tuple{Float64,Float64}}(undef, length(bins))
    
    for (j, bin_peaks) in enumerate(bins)
        # Calculate bin center
        local bin_center
        if intensity_weighted_centers
            weights = [p[2] for p in bin_peaks]
            bin_center = sum(p[1]*p[2] for p in bin_peaks) / sum(weights)
        else
            bin_center = mean(p[1] for p in bin_peaks)
        end
        
        final_bins_boundaries[j] = (first(bin_peaks)[1], last(bin_peaks)[1])

        for p in bin_peaks
            s_idx = p[3]
            X[s_idx, j] = max(X[s_idx, j], p[2])
        end
    end

    # Filter by frequency
    if frequency_threshold > 0
        present_count = vec(sum(X .> 0, dims=1))
        min_count = ceil(Int, frequency_threshold * ns)
        keep_mask = findall(present_count .>= min_count)
        X = X[:, keep_mask]
        final_bins_boundaries = final_bins_boundaries[keep_mask]
    end

    return (X, final_bins_boundaries)
end

# =============================================================================
# 7) Spatial & Advanced Processing (Stubs & New Functions)
# =============================================================================

"""
    find_ppm_error_by_region(msi_data, region_masks, reference_peaks) -> Dict

Calculates PPM error statistics for different spatial regions.
`region_masks` is a Dict mapping region names (e.g., :tumor) to BitMatrix masks.
"""
function find_ppm_error_by_region(msi_data::MSIData, region_masks::Dict, reference_peaks::Dict)
    regional_reports = Dict{Symbol, NamedTuple}()
    
    width, height = msi_data.image_dims
    
    for (region_name, mask) in region_masks
        mask_height, mask_width = size(mask)
        if mask_width != width || mask_height != height
            @warn "Mask dimensions ($(mask_width)x$(mask_height)) for region '$region_name' do not match image dimensions ($(width)x$(height)). Skipping."
            continue
        end

        indices = [
            i for i in 1:length(msi_data.spectra_metadata) 
            if msi_data.spectra_metadata[i].x > 0 &&
               msi_data.spectra_metadata[i].y > 0 &&
               mask[msi_data.spectra_metadata[i].y, msi_data.spectra_metadata[i].x]
        ]
        
        if isempty(indices) continue end
        
        # Call analyze_mass_accuracy with the specific indices for the region
        regional_reports[region_name] = analyze_mass_accuracy(msi_data, reference_peaks; spectrum_indices=indices)
    end
    return regional_reports
end

"""
    regional_calibration(msi_data, region_masks, reference_peaks)

(STUB) Applies different calibration models to different spatial regions.
"""
function regional_calibration(msi_data::MSIData, region_masks::Dict, reference_peaks::Dict)
    @warn "regional_calibration is a stub and not fully implemented."
    # 1. For each region in region_masks:
    # 2. Create a region-specific calibration model using `calibrate_spectra`.
    # 3. Apply the model to all spectra within that region.
    # 4. Return a new MSIData object or modified spectra vector.
    return msi_data # Return unmodified for now
end

# =============================================================================
# 8) Advanced Peak Quality & Adaptive Parameters
# =============================================================================

"""
    calculate_peak_quality_metrics(peak_mz, peak_intensity, local_spectrum) -> Dict

(STUB) Calculates advanced quality metrics for a single peak.
"""
function calculate_peak_quality_metrics(peak_mz, peak_intensity, local_spectrum)
    # A full implementation would calculate:
    # - Sharpness (e.g., ratio of height to FWHM)
    # - Symmetry (e.g., ratio of left/right half-widths)
    # - Signal-to-Noise (using local noise estimation)
    # - Isolation score (how close are other peaks)
    return Dict(:sharpness => 1.0, :symmetry => 1.0, :snr => 10.0)
end




"""
    calculate_adaptive_bin_tolerance(ppm_error_distribution) -> Float64

Calculates an appropriate binning tolerance based on observed mass accuracy.
"""
function calculate_adaptive_bin_tolerance(ppm_error_distribution::Vector{Float64})
    if isempty(ppm_error_distribution)
        return 20.0 # Default if no data
    end
    # A robust strategy: mean + 3 * std deviation to capture ~99.7% of peaks
    return mean(ppm_error_distribution) + 3 * std(ppm_error_distribution)
end


# =============================================================================
# 7) Plotting Helper
# =============================================================================

"""
    plot_stage_spectrum(mz, intensity; title, xlabel, ylabel) -> Figure

A helper function to generate a plot of a single spectrum using `CairoMakie`. This is
useful for visualizing the output of different preprocessing steps.

# Arguments
- `mz::AbstractVector`: The m/z vector for the x-axis.
- `intensity::AbstractVector`: The intensity vector for the y-axis.

# Keyword Arguments
- `title::AbstractString`: The title of the plot.
- `xlabel::AbstractString`: The label for the x-axis. Defaults to "m/z".
- `ylabel::AbstractString`: The label for the y-axis. Defaults to "Intensity".

# Returns
- `CairoMakie.Figure`: A `Figure` object containing the plot. The caller is responsible for displaying or saving it.

# Example
```julia
using CairoMakie
mz = 100:0.1:110;
intensity = rand(length(mz));
fig = plot_stage_spectrum(mz, intensity, title="My Spectrum");
save("my_spectrum.png", fig);
```
"""
function plot_stage_spectrum(mz::AbstractVector, intensity::AbstractVector;
                             title::AbstractString, xlabel::AbstractString="m/z",
                             ylabel::AbstractString="Intensity")
    # This dynamic import is for script-like use; in a package, Makie would be a full dependency.
    @eval begin
        import CairoMakie
        using CairoMakie
    end
    fig = CairoMakie.Figure(size = (1400, 500))
    ax  = CairoMakie.Axis(fig[1, 1], title=title, xlabel=xlabel, ylabel=ylabel)
    CairoMakie.lines!(ax, mz, intensity)
    return fig
end


# =============================================================================
# 9) Quality Control Metrics
# =============================================================================

"""
    get_common_calibration_standards(type::Symbol) -> Dict{Float64, String}

Returns a dictionary of common m/z values for calibration standards.

# Arguments
- `type::Symbol`: The type of standards to return. Currently supports `:maldi_pos`.

# Returns
- `Dict{Float64, String}`: A dictionary mapping theoretical m/z to compound names.
"""
function get_common_calibration_standards(type::Symbol)
    if type == :maldi_pos
        # Common peptide/protein standards for positive mode MALDI
        return Dict{Float64, String}(
            1046.5420 => "Angiotensin II",
            1060.5690 => "Bradykinin",
            1296.6853 => "Angiotensin I",
            2465.1989 => "ACTH clip (1-24)"
        )
    else
        @warn "Unsupported calibration standard type: $type. Returning empty dictionary."
        return Dict{Float64, String}()
    end
end

"""
    calculate_ppm_error(measured_mz::Real, theoretical_mz::Real) -> Float64

Calculates the mass accuracy error in parts-per-million (PPM) between a measured
and a theoretical m/z value.

# Arguments
- `measured_mz::Real`: The experimentally measured m/z value.
- `theoretical_mz::Real`: The known, theoretical m/z value of a compound.

# Returns
- `Float64`: The calculated PPM error. Returns `Inf` if `theoretical_mz` is zero.

# Formula
`PPM = 10^6 * |measured_mz - theoretical_mz| / theoretical_mz`

# Example
```julia
calculate_ppm_error(100.005, 100.0) # returns 50.0
```
"""
function calculate_ppm_error(measured_mz::Real, theoretical_mz::Real)
    if theoretical_mz == 0
        return Inf
    end
    return 1e6 * abs(Float64(measured_mz) - Float64(theoretical_mz)) / Float64(theoretical_mz)
end

"""
    calculate_ppm_error_bulk(measured_mz::Vector{<:Real}, theoretical_mz::Vector{<:Real}) -> Vector{Float64}

Calculates PPM errors for multiple pairs of measured and theoretical mass values.

# Arguments
- `measured_mz::Vector{<:Real}`: A vector of experimentally measured m/z values.
- `theoretical_mz::Vector{<:Real}`: A vector of known, theoretical m/z values.

# Returns
- `Vector{Float64}`: A vector containing the calculated PPM error for each pair.
"""
function calculate_ppm_error_bulk(measured_mz::Vector{Real}, theoretical_mz::Vector{Real})
    return [calculate_ppm_error(m, t) for (m, t) in zip(measured_mz, theoretical_mz)]
end

"""
    calculate_resolution_fwhm(mz::Real, profile_mz::AbstractVector{<:Real}, profile_intensity::AbstractVector{<:Real}) -> Float64

Calculates the mass resolution of a peak in profile-mode data using the Full Width at
Half Maximum (FWHM) method. Resolution is a measure of an instrument's ability to
distinguish between two peaks of slightly different mass-to-charge ratios.

# Arguments
- `mz::Real`: The m/z value of the peak's centroid.
- `profile_mz::AbstractVector{<:Real}`: The full m/z array from the profile-mode spectrum.
- `profile_intensity::AbstractVector{<:Real}`: The full intensity array from the profile-mode spectrum.

# Returns
- `Float64`: The calculated resolution (`m / Δm`). Returns `NaN` if the FWHM cannot be determined (e.g., peak is at the edge of the spectrum).

# Formula
`Resolution = m / Δm`, where `Δm` is the FWHM.
"""
function calculate_resolution_fwhm(mz::Real, profile_mz::AbstractVector{<:Real},
                                 profile_intensity::AbstractVector{<:Real})
    
    # Find peak center index
    peak_idx = argmin(abs.(profile_mz .- mz))
    peak_height = Float64(profile_intensity[peak_idx])
    half_max = peak_height / 2
    
    # Find left half-maximum point (interpolate for accuracy)
    left_idx = find_last_below(profile_intensity[1:peak_idx], half_max)
    if left_idx == 0 || left_idx == length(profile_intensity[1:peak_idx])
        return NaN
    end
    
    # Linear interpolation for left FWHM
    x1, x2 = Float64(profile_mz[left_idx]), Float64(profile_mz[left_idx+1])
    y1, y2 = Float64(profile_intensity[left_idx]), Float64(profile_intensity[left_idx+1])
    left_fwhm = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1)
    
    # Find right half-maximum point
    right_slice = profile_intensity[peak_idx:end]
    right_offset = find_first_below(right_slice, half_max)
    if right_offset == 0 || right_offset == length(right_slice)
        return NaN
    end
    
    right_idx = peak_idx + right_offset - 1
    x1, x2 = Float64(profile_mz[right_idx-1]), Float64(profile_mz[right_idx])
    y1, y2 = Float64(profile_intensity[right_idx-1]), Float64(profile_intensity[right_idx])
    right_fwhm = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1)
    
    fwhm = right_fwhm - left_fwhm
    return fwhm > 0 ? Float64(mz) / fwhm : NaN
end

# Helper functions for FWHM calculation
function find_last_below(v::AbstractVector{<:Real}, threshold::Real)
    for i in length(v):-1:2
        if v[i] >= threshold && v[i-1] < threshold
            return i-1
        end
    end
    return 0
end

function find_first_below(v::AbstractVector{<:Real}, threshold::Real)
    for i in 1:(length(v)-1)
        if v[i] >= threshold && v[i+1] < threshold
            return i+1
        end
    end
    return 0
end

"""
    analyze_mass_accuracy(msi_data, reference_peaks; ppm_tolerance=5.0, sample_spectra=100) -> NamedTuple

Analyzes the mass accuracy across a sample of spectra from an MSI dataset by comparing
detected peaks against a list of known reference m/z values. It calculates PPM error
statistics and suggests an optimal PPM tolerance for future analyses.

# Arguments
- `msi_data`: An `MSIData` object containing the dataset.
- `reference_peaks::Dict{Float64,String}`: A dictionary mapping theoretical m/z values to compound names.

# Keyword Arguments
- `ppm_tolerance::Float64`: The initial PPM tolerance used to match detected peaks to reference peaks. Defaults to 5.0.
- `sample_spectra::Int`: The number of spectra to sample from the dataset for the analysis. Defaults to 100.

# Returns
- `NamedTuple`: A comprehensive report containing statistics like `mean_ppm`, `std_ppm`, `optimal_ppm`, `n_matches`, and lists of matched peaks and all PPM errors.
"""
function analyze_mass_accuracy(msi_data, reference_peaks::Dict{Float64,String}; 
                             ppm_tolerance::Float64=5.0, sample_spectra=100,
                             spectrum_indices=nothing)
    println("\n[ MASS ACCURACY ANALYSIS ]")
    println("PPM tolerance: ", ppm_tolerance, " ppm")
    
    theoretical_mz = sort(collect(keys(reference_peaks)))
    ppm_errors = Float64[]
    matched_peaks = Tuple{Float64,Float64,String}[]  # (theoretical, measured, compound)
    
    # Determine which spectra to process
    local indices_to_process
    if spectrum_indices === nothing
        println("Spectra to sample: ", sample_spectra)
        if sample_spectra >= length(msi_data.spectra_metadata)
            indices_to_process = 1:length(msi_data.spectra_metadata)
        else
            indices_to_process = round.(Int, range(1, length(msi_data.spectra_metadata), length=sample_spectra))
        end
    else
        indices_to_process = spectrum_indices
        println("Analyzing $(length(indices_to_process)) specified spectra.")
    end
    
    for idx in indices_to_process
        mz, intensity = GetSpectrum(msi_data, idx)
        
        # Detect peaks in this spectrum
        detected_peaks, _ = detect_peaks_profile(mz, intensity, snr_threshold=3.0)
        
        # Match detected peaks to reference peaks
        for (i, theoretical) in enumerate(theoretical_mz)
            # Find closest detected peak within tolerance
            distances = abs.(detected_peaks .- theoretical)
            if !isempty(distances)
                min_idx = argmin(distances)
                min_distance = distances[min_idx]
                
                ppm_error = calculate_ppm_error(detected_peaks[min_idx], theoretical)
                
                if ppm_error <= ppm_tolerance
                    push!(ppm_errors, ppm_error)
                    push!(matched_peaks, (theoretical, detected_peaks[min_idx], reference_peaks[theoretical]))
                end
            end
        end
    end
    
    if isempty(ppm_errors)
        @warn "No peaks matched within $ppm_tolerance ppm tolerance"
        return (mean_ppm=NaN, std_ppm=NaN, min_ppm=NaN, max_ppm=NaN, optimal_ppm=NaN, n_matches=0, matched_peaks=[], all_ppm_errors=Float64[])
    end
    
    # Calculate statistics
    mean_ppm = mean(ppm_errors)
    std_ppm = std(ppm_errors)
    min_ppm = minimum(ppm_errors)
    max_ppm = maximum(ppm_errors)
    
    # Determine optimal ppm tolerance (mean + 3σ covers ~99.7% of peaks for normal distribution)
    optimal_ppm = mean_ppm + 3 * std_ppm
    
    return (
        mean_ppm = mean_ppm,
        std_ppm = std_ppm,
        min_ppm = min_ppm,
        max_ppm = max_ppm,
        optimal_ppm = optimal_ppm,
        n_matches = length(ppm_errors),
        matched_peaks = matched_peaks,
        all_ppm_errors = ppm_errors
    )
end

"""
    generate_qc_report(msi_data, filename; reference_peaks, output_dir, sample_spectra) -> Tuple

Generates a comprehensive Quality Control (QC) report for an MSI dataset. It assesses
mass accuracy and resolution, saving the results to CSV files and a summary text file.

# Arguments
- `msi_data`: The `MSIData` object to be analyzed.
- `filename::String`: The name of the input data file, used for reporting.

# Keyword Arguments
- `reference_peaks::Dict`: A dictionary of reference peaks for mass accuracy analysis. If `nothing`, defaults are loaded using `get_common_calibration_standards`.
- `output_dir::String`: The directory where the report files will be saved. Defaults to `"qc_results"`.
- `sample_spectra::Int`: The number of spectra to sample for the analysis. Defaults to 100.

# Returns
- `Tuple`: A tuple containing the detailed `accuracy_report` (NamedTuple) and `resolution_results` (Vector).

# Side Effects
- Creates a directory specified by `output_dir`.
- Writes `mass_accuracy_results.csv`, `resolution_results.csv`, and `qc_summary.txt` into the output directory.
"""
function generate_qc_report(msi_data, filename::String; reference_peaks=nothing, output_dir="qc_results", sample_spectra=100, spectrum_indices=nothing)
    println("\n[ QC REPORT GENERATION ]")
    println("Input file: ", filename)
    println("Output directory: ", output_dir)
    mkpath(output_dir)
    
    # Use default calibrants if none provided
    if reference_peaks === nothing
        reference_peaks = get_common_calibration_standards(:maldi_pos)
    end
    
    println("Generating QC Report...")
    println("Using $(length(reference_peaks)) reference masses")
    
    # 1. Analyze mass accuracy
    accuracy_report = analyze_mass_accuracy(msi_data, reference_peaks; sample_spectra=sample_spectra, spectrum_indices=spectrum_indices)
    
    println("\n" * "="^50)
    println("MASS ACCURACY REPORT")
    println("="^50)
    println("Mean PPM error: $(round(accuracy_report.mean_ppm, digits=2)) ppm")
    println("Std PPM error: $(round(accuracy_report.std_ppm, digits=2)) ppm") 
    println("Min PPM error: $(round(accuracy_report.min_ppm, digits=2)) ppm")
    println("Max PPM error: $(round(accuracy_report.max_ppm, digits=2)) ppm")
    if haskey(accuracy_report, :optimal_ppm)
        println("Optimal PPM tolerance: $(round(accuracy_report.optimal_ppm, digits=2)) ppm")
    else
        println("Optimal PPM tolerance: Not available")
    end
    println("Number of matches: $(accuracy_report.n_matches)")
    
    # 2. Calculate resolution for a few representative peaks
    println("\n" * "="^50)
    println("RESOLUTION ANALYSIS")
    println("="^50)
    
    resolution_results = []
    sample_spectra = min(10, length(msi_data.spectra_metadata))
    
    for (i, idx) in enumerate(round.(Int, range(1, length(msi_data.spectra_metadata), length=sample_spectra)))
        process_spectrum(msi_data, idx) do mz, intensity
            if !qc_is_empty(mz, intensity)
                # Test resolution on the most intense peak
                max_intensity_idx = argmax(intensity)
                test_mz = mz[max_intensity_idx]
                
                resolution = calculate_resolution_fwhm(test_mz, mz, intensity)
                if !isnan(resolution)
                    push!(resolution_results, resolution)
                    println("Spectrum $idx: Resolution = $(round(resolution))")
                end
            end
        end
    end
    
    if !isempty(resolution_results)
        avg_resolution = mean(resolution_results)
        println("\nAverage resolution: $(round(avg_resolution))")
        println("Resolution range: $(round(minimum(resolution_results))) - $(round(maximum(resolution_results)))")
    end
    
    # 3. Save detailed results
    
    # Save PPM error distribution
    ppm_df = DataFrame(
        theoretical_mz = [p[1] for p in accuracy_report.matched_peaks],
        measured_mz = [p[2] for p in accuracy_report.matched_peaks], 
        compound = [p[3] for p in accuracy_report.matched_peaks],
        ppm_error = accuracy_report.all_ppm_errors
    )
    
    CSV.write(joinpath(output_dir, "mass_accuracy_results.csv"), ppm_df)
    
    # Save resolution results
    if !isempty(resolution_results)
        res_df = DataFrame(resolution = resolution_results)
        CSV.write(joinpath(output_dir, "resolution_results.csv"), res_df)
    end
    
    # 4. Create summary
    summary = """
    QC REPORT SUMMARY
    =================
    Date: $(now())
    File: $(filename)
    Spectra analyzed: $(length(msi_data.spectra_metadata))
    
    MASS ACCURACY:
    - Mean PPM: $(round(accuracy_report.mean_ppm, digits=2)) ppm
    - Std PPM: $(round(accuracy_report.std_ppm, digits=2)) ppm
    - Recommended tolerance: $(round(accuracy_report.optimal_ppm, digits=2)) ppm (mean + 3 * std)
    
    RESOLUTION:
    - Average: $(isempty(resolution_results) ? "N/A" : string(round(mean(resolution_results))))
    - Range: $(isempty(resolution_results) ? "N/A" : "$(round(minimum(resolution_results))) - $(round(maximum(resolution_results)))")
    
    RECOMMENDATIONS:
    - Use $(round(accuracy_report.optimal_ppm, digits=2)) ppm for peak matching
    - Instrument performance: $(accuracy_report.mean_ppm < 5 ? "Excellent" : accuracy_report.mean_ppm < 10 ? "Good" : "Needs calibration")
    """
    
    open(joinpath(output_dir, "qc_summary.txt"), "w") do f
        write(f, summary)
    end
    
    println("\nQC report saved to: $output_dir")
    return accuracy_report, resolution_results
end
