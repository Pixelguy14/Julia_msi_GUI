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

# =============================================================================
# Data Structures
# =============================================================================

"""
    FeatureMatrix

Structure to hold the final feature matrix, including m/z bin boundaries and sample indices.

# Fields
- `matrix`: The numerical matrix where rows are samples and columns are features (m/z bins).
- `mz_bins`: A vector of tuples `(low_mz, high_mz)` for each feature column.
- `sample_ids`: A vector of indices corresponding to the original spectra.
"""
struct FeatureMatrix
    matrix::Array{Float64,2}                 # samples × features
    mz_bins::Vector{Tuple{Float64,Float64}}  # [(low, high), ...]
    sample_ids::Vector{Int}
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

Checks that the m/z axis is monotonically non-decreasing, as expected in a profile spectrum.
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
Supported methods: `:sqrt` (default) and `:log1p`.
"""
function transform_intensity(intensity::AbstractVector{<:Real}; method::Symbol=:sqrt)
    if method === :sqrt
        return sqrt.(max.(zero(eltype(intensity)), intensity))
    elseif method === :log1p
        return log1p.(max.(zero(eltype(intensity)), intensity))
    else
        return collect(float.(intensity))
    end
end

"""
    smooth_spectrum(y; window=21, order=2) -> Vector

Applies a Savitzky–Golay filter if `SavitzkyGolay.jl` is available.
Otherwise, falls back to a simple (non-phase-correct) moving average.
"""
function smooth_spectrum(y::AbstractVector{<:Real}; window::Int=9, order::Int=2)
    win = isodd(window) ? window : window + 1
    res = SavitzkyGolay.savitzky_golay(collect(float.(y)), win, order)
    return res.y
end

# =============================================================================
# 2) Baseline Correction
# =============================================================================

"""
    snip_baseline(y, iterations=100) -> Vector

Estimates the baseline using a simple 1D SNIP (Statistics-sensitive Non-linear
Iterative Peak-clipping) algorithm. `iterations` controls the aggressiveness.
"""
function snip_baseline(y::AbstractVector{<:Real}, iterations::Int=100)
    n = length(y)
    b = collect(float.(y))  # work copy
    buf = similar(b)
    for k in 1:iterations
        copyto!(buf, b)
        @inbounds for i in 2:n-1
            buf[i] = min(b[i], 0.5 * (b[i-1] + b[i+1]))
        end
        # Handle endpoints
        buf[1] = min(b[1], b[2])
        buf[end] = min(b[end], b[end-1])
        b, buf = buf, b # Swap buffers
    end
    return b
end

# =============================================================================
# 3) Intensity Normalization
# =============================================================================

"""
    tic_normalize(y) -> Vector

Normalizes intensities to the Total Ion Current (TIC). If the sum is zero, returns a copy.
"""
function tic_normalize(y::AbstractVector{<:Real})
    s = sum(y)
    return s <= 0 ? collect(float.(y)) : collect(float.(y)) ./ s
end

"""
    pqn_normalize(M) -> Matrix

Performs Probabilistic Quotient Normalization on a matrix `M` where columns are spectra.
"""
function pqn_normalize(M::AbstractMatrix{<:Real})
    M_float = collect(float.(M))
    # Calculate reference spectrum (median across all spectra)
    ref = mapslices(median, M_float; dims=2)[:,1]
    
    # Calculate quotients for each spectrum relative to the reference
    Q = similar(M_float)
    @inbounds for j in axes(M_float, 2)
        Q[:, j] = M_float[:, j] ./ (ref .+ eps(eltype(M_float)))
    end
    
    # Find the median quotient for each spectrum (scaling factor)
    s = [median( @view Q[:, j]) for j in axes(Q, 2)]
    
    # Normalize the original matrix
    @inbounds for j in axes(M_float, 2)
        M_float[:, j] ./= (s[j] + eps(eltype(M_float)))
    end
    return M_float
end

# =============================================================================
# 4) Peak Detection (for Profile Data)
# =============================================================================

"""
    detect_peaks_profile(mz, y; half_window=10, snr_threshold=2.0)

Detects local maxima with a signal-to-noise threshold (using MAD for noise estimation).
Assumes profile-mode data and a monotonic m/z axis.
"""
function detect_peaks_profile(mz::AbstractVector{<:Real},
                              y::AbstractVector{<:Real};
                              half_window::Int=10,
                              snr_threshold::Float64=2.0)
    n = length(y)
    n < 3 && return (Float64[], Float64[])
    
    # Noise estimation using Median Absolute Deviation (robust to peaks)
    noise_level = mad(y, normalize=true) + eps(Float64)

    # Smooth the spectrum to make peak detection more robust
    ys = smooth_spectrum(y; window=max(5, 2*half_window+1), order=2)

    peak_idx = Int[]
    @inbounds for i in 2:n-1
        left = max(1, i - half_window)
        right = min(n, i + half_window)
        
        local_max = ys[i]
        # A point is a peak if it's the maximum in its neighborhood and above the SNR threshold
        if local_max >= maximum( @view ys[left:right]) && (local_max > snr_threshold * noise_level)
            # Ensure we only record one point for flat-topped peaks
            if isempty(peak_idx) || (i - last(peak_idx) > half_window)
                push!(peak_idx, i)
            end
        end
    end

    # Return original intensities at peak locations
    pk_mz  = [float(mz[i]) for i in peak_idx]
    pk_int = [float(y[i]) for i in peak_idx]
    
    return (pk_mz, pk_int)
end

"""
    detect_peaks_centroid(mz, y; intensity_threshold=0.0)

Filters centroided data based on a minimum intensity threshold.
"""
function detect_peaks_centroid(mz::AbstractVector{<:Real},
                               y::AbstractVector{<:Real};
                               intensity_threshold::Float64=0.0)
    
    keep_indices = findall(y .>= intensity_threshold)
    
    return (mz[keep_indices], y[keep_indices])
end

# =============================================================================
# 5) Peak Alignment
# =============================================================================

"""
    align_peaks_lowess(ref_mz, tgt_mz; tolerance=0.002) -> warp::Function

Generates a warping function `warp(x)` to map target m/z values to reference m/z values.
Uses a lightweight LOWESS-like approach with linear interpolation.
"""
function align_peaks_lowess(ref_mz::Vector{<:Real},
                            tgt_mz::Vector{<:Real};
                            tolerance::Float64=0.002)
    # Efficiently match peaks between sorted lists
    pairs = Tuple{Float64,Float64}[]  # (target_mz, reference_mz)
    i = 1; j = 1
    while i <= length(tgt_mz) && j <= length(ref_mz)
        dt = tgt_mz[i] - ref_mz[j]
        if abs(dt) <= tolerance
            push!(pairs, (float(tgt_mz[i]), float(ref_mz[j])))
            i += 1; j += 1
        elseif dt < 0
            i += 1
        else
            j += 1
        end
    end

    if length(pairs) < 3
        @warn "Too few matching peaks for alignment. Returning identity function."
        return x -> float.(x)
    end

    t = [p[1] for p in pairs]
    r = [p[2] for p in pairs]

    # Lightly smooth the mapping to reduce noise
    t_s = smooth_spectrum(t; window=5, order=2)
    r_s = smooth_spectrum(r; window=5, order=2)

    # Return a function that performs linear interpolation for warping
    function warp(x::AbstractVector{<:Real})
        out = similar(collect(float.(x)))
        for (k, xv) in enumerate(x)
            if xv <= t_s[1]
                # Linear extrapolation at the start
                m = (r_s[2]-r_s[1]) / (t_s[2]-t_s[1] + eps())
                out[k] = r_s[1] + m*(xv - t_s[1])
            elseif xv >= t_s[end]
                # Linear extrapolation at the end
                m = (r_s[end]-r_s[end-1]) / (t_s[end]-t_s[end-1] + eps())
                out[k] = r_s[end-1] + m*(xv - t_s[end-1])
            else
                # Linear interpolation for points in the middle
                lo = searchsortedlast(t_s, xv)
                hi = lo + 1
                α = (xv - t_s[lo]) / (t_s[hi] - t_s[lo] + eps())
                out[k] = (1-α)*r_s[lo] + α*r_s[hi]
            end
        end
        return out
    end

    return warp
end

# =============================================================================
# 6) Peak Binning & Feature Matrix Generation
# =============================================================================

"""
    _find_bin_index(x, bins) -> Int

Efficiently finds the index of the bin `(low, high)` that contains `x` using binary search.
Returns 0 if not found.
"""
function _find_bin_index(x::Float64, bins::Vector{Tuple{Float64,Float64}})
    lo, hi = 1, length(bins)
    while lo <= hi
        mid = (lo + hi) >>> 1
        b = bins[mid]
        if x < b[1]
            hi = mid - 1
        elseif x > b[2]
            lo = mid + 1
        else
            return mid
        end
    end
    return 0
end

"""
    bin_peaks(all_pk_mz, all_pk_int, tolerance; frequency_threshold=0.25)

Groups peaks from all spectra into consensus m/z bins and creates a feature matrix.
Filters out features that do not appear in a minimum fraction of spectra.
"""
function bin_peaks(all_pk_mz::Vector{<:AbstractVector{<:Real}},
                   all_pk_int::Vector{<:AbstractVector{<:Real}},
                   tolerance::Float64; frequency_threshold::Float64=0.25)
    ns = length(all_pk_mz)
    ns == 0 && return (zeros(0,0), Tuple{Float64,Float64}[])

    # 1) Collect all unique peak m/z values and sort them
    flat_mz = Float64[]
    for v in all_pk_mz
        append!(flat_mz, float.(v))
    end
    sort!(flat_mz)
    isempty(flat_mz) && return (zeros(ns, 0), Tuple{Float64,Float64}[])

    # 2) Create contiguous m/z bins based on tolerance
    bins = Tuple{Float64,Float64}[]
    cur_lo = flat_mz[1]
    cur_hi = flat_mz[1]
    for x in @view flat_mz[2:end]
        if x - cur_hi <= tolerance
            cur_hi = x # Extend the current bin
        else
            push!(bins, (cur_lo, cur_hi)) # Finalize old bin
            cur_lo = x; cur_hi = x      # Start a new one
        end
    end
    push!(bins, (cur_lo, cur_hi))

    # 3) Create the feature matrix (samples x features) using max intensity per bin
    X = zeros(Float64, ns, length(bins))
    for i in 1:ns
        for (mzv, iv) in zip(all_pk_mz[i], all_pk_int[i])
            bidx = _find_bin_index(float(mzv), bins)
            if bidx > 0
                X[i, bidx] = max(X[i, bidx], float(iv))
            end
        end
    end

    # 4) Filter features by minimum frequency
    if frequency_threshold > 0
        present_count = vec(sum(X .> 0, dims=1))
        min_count = ceil(Int, frequency_threshold * ns)
        keep_mask = findall(present_count .>= min_count)
        X = X[:, keep_mask]
        bins = bins[keep_mask]
    end

    return (X, bins)
end

# =============================================================================
# 7) Plotting Helper
# =============================================================================

"""
    plot_stage_spectrum(mz, intensity; title, ...)

Returns a `CairoMakie.Figure` for a single spectrum trace. The caller is responsible for saving.
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
# 8) Pipeline Orchestrator
# =============================================================================

"""
    run_preprocessing_pipeline(spectra; steps, params, on_stage)

Executes a flexible preprocessing pipeline on a vector of spectra.

# Arguments
- `spectra`: A vector of `(mz, intensity)` tuples.
- `steps`: A vector of symbols defining the pipeline order (e.g., `[:qc, :smooth, :baseline, :peaks, :bin]`).
- `params`: A dictionary of parameters for each step.
- `on_stage`: An optional callback function `on_stage(stage_symbol; idx, mz, intensity)` executed after each step for logging or visualization.

# Returns
- A `FeatureMatrix` if `:bin` is in the steps, otherwise the vector of processed spectra.
"""
function run_preprocessing_pipeline(spectra::Vector;
                                    steps::Vector{Symbol},
                                    params::Dict=Dict(),
                                    on_stage::Function=(;kwargs...)->nothing)

    processed = deepcopy(spectra)  # Don't mutate the original input
    reference_peaks = nothing      # For alignment

    _emit(stage::Symbol, idx::Int, mz, y) = on_stage(stage; idx=idx, mz=mz, intensity=y)

    for step in steps
        @info "Running step: $step"

        if step === :qc
            for (i, (mz, y)) in enumerate(processed)
                (isempty(mz) || isempty(y)) && continue
                _emit(:qc_raw, i, mz, y)
                qc_is_empty(mz, y) && @warn "Spectrum at index $i is empty."
                !qc_is_regular(mz) && @warn "m/z axis at index $i is not monotonic."
            end

        elseif step === :transform
            meth = get(params, :transform_method, :sqrt)
            for i in eachindex(processed)
                mz, y = processed[i]
                y_new = transform_intensity(y; method=meth)
                processed[i] = (mz, y_new)
                _emit(:transform, i, mz, y_new)
            end

        elseif step === :smooth
            win = get(params, :sg_window, 21)
            ord = get(params, :sg_order, 2)
            for i in eachindex(processed)
                mz, y = processed[i]
                y_smooth = smooth_spectrum(y; window=win, order=ord)
                processed[i] = (mz, y_smooth)
                _emit(:smooth, i, mz, y_smooth)
            end

        elseif step === :baseline
            iters = get(params, :snip_iterations, 100)
            for i in eachindex(processed)
                mz, y = processed[i]
                baseline = snip_baseline(y, iters)
                y_corrected = max.(0.0, y .- baseline)
                processed[i] = (mz, y_corrected)
                _emit(:baseline, i, mz, y_corrected)
            end

        elseif step === :normalize
            mode = get(params, :normalize_method, :tic)
            if mode === :tic
                for i in eachindex(processed)
                    mz, y = processed[i]
                    y_norm = tic_normalize(y)
                    processed[i] = (mz, y_norm)
                    _emit(:normalize, i, mz, y_norm)
                end
            elseif mode === :pqn
                # Note: PQN assumes spectra are on a common m/z grid.
                matrix = hcat([float.(p[2]) for p in processed]...)
                matrix_norm = pqn_normalize(matrix)
                for i in eachindex(processed)
                    mz, _ = processed[i]
                    processed[i] = (mz, view(matrix_norm, :, i))
                    _emit(:normalize, i, mz, processed[i][2])
                end
            end

        elseif step === :peaks
            peak_results = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, length(processed))
            hw = get(params, :peak_half_window, 10)
            snr = get(params, :peak_snr, 2.0)
            for (i, (mz, y)) in enumerate(processed)
                pk_mz, pk_int = detect_peaks_profile(mz, y; half_window=hw, snr_threshold=snr)
                peak_results[i] = (pk_mz, pk_int)
                # Emit with original mz axis but maybe stem plot of peaks?
                _emit(:peaks, i, pk_mz, pk_int)
            end
            processed = peak_results
            # Set reference for alignment
            reference_peaks = isempty(processed) ? nothing : processed[1][1]

        elseif step === :align
            reference_peaks === nothing && (@error "Alignment requires a :peaks step first."; continue)
            tol = get(params, :align_tolerance, 0.002)
            for i in 2:length(processed)
                tgt_peaks, intens = processed[i]
                warp_func = align_peaks_lowess(reference_peaks, tgt_peaks; tolerance=tol)
                processed[i] = (warp_func(tgt_peaks), intens)
                _emit(:align, i, processed[i][1], processed[i][2])
            end

        elseif step === :bin
            all_pks  = [s[1] for s in processed]
            all_ints = [s[2] for s in processed]
            tol  = get(params, :bin_tolerance, 0.002)
            freq = get(params, :bin_min_frequency, 0.25)
            mat, mz_bins = bin_peaks(all_pks, all_ints, tol; frequency_threshold=freq)
            return FeatureMatrix(mat, mz_bins, collect(1:length(processed)))
        end
    end

    @warn "Pipeline finished without a :bin step; returning processed spectra."
    return processed
end

"""
    run_preprocessing_pipeline(msi_data::MSIData, indices::Vector{Int}; steps, params, on_stage)

Executes a flexible preprocessing pipeline on a subset of spectra from an MSIData object,
with mode-aware logic for centroid and profile data.
"""
function run_preprocessing_pipeline(msi_data::MSIData, indices::Vector{Int};
                                    steps::Vector{Symbol},
                                    params::Dict=Dict(),
                                    on_stage::Function=(;kwargs...)->nothing)

    # This version of the pipeline is mode-aware.
    # It processes spectra directly from the MSIData object.

    processed_spectra = Vector{Tuple}(undef, length(indices))

    # First, load all spectra and apply initial steps that run on individual spectra
    for (i, spec_idx) in enumerate(indices)
        mz, intensity = GetSpectrum(msi_data, spec_idx)
        mode = msi_data.spectra_metadata[spec_idx].mode

        on_stage(:qc_raw; idx=spec_idx, mz=mz, intensity=intensity)

        for step in steps
            if step === :transform
                meth = get(params, :transform_method, :sqrt)
                intensity = transform_intensity(intensity; method=meth)
                on_stage(:transform; idx=spec_idx, mz=mz, intensity=intensity)
            elseif step === :smooth && mode == PROFILE
                win = get(params, :sg_window, 21)
                ord = get(params, :sg_order, 2)
                intensity = smooth_spectrum(intensity; window=win, order=ord)
                on_stage(:smooth; idx=spec_idx, mz=mz, intensity=intensity)
            elseif step === :baseline && mode == PROFILE
                iters = get(params, :snip_iterations, 100)
                baseline = snip_baseline(intensity, iters)
                intensity = max.(0.0, intensity .- baseline)
                on_stage(:baseline; idx=spec_idx, mz=mz, intensity=intensity)
            elseif step === :normalize
                norm_mode = get(params, :normalize_method, :tic)
                if norm_mode === :tic
                    intensity = tic_normalize(intensity)
                    on_stage(:normalize; idx=spec_idx, mz=mz, intensity=intensity)
                end
            end
        end
        processed_spectra[i] = (mz, intensity, mode) # Store mode for peak detection
    end

    # Now, handle steps that require all spectra (like PQN) or are the final steps
    final_result = nothing
    for step in steps
        if step === :normalize && get(params, :normalize_method, :tic) === :pqn
            # Note: PQN assumes spectra are on a common m/z grid.
            matrix = hcat([float.(p[2]) for p in processed_spectra]...)
            matrix_norm = pqn_normalize(matrix)
            for i in eachindex(processed_spectra)
                mz, _, mode = processed_spectra[i]
                processed_spectra[i] = (mz, view(matrix_norm, :, i), mode)
                on_stage(:normalize; idx=indices[i], mz=mz, intensity=processed_spectra[i][2])
            end
        elseif step === :peaks
            peak_results = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, length(processed_spectra))
            hw = get(params, :peak_half_window, 10)
            snr = get(params, :peak_snr, 2.0)
            intensity_thresh = get(params, :peak_intensity_threshold, 0.0)

            for (i, (mz, y, mode)) in enumerate(processed_spectra)
                spec_idx = indices[i]
                if mode == PROFILE
                    pk_mz, pk_int = detect_peaks_profile(mz, y; half_window=hw, snr_threshold=snr)
                else # CENTROID
                    pk_mz, pk_int = detect_peaks_centroid(mz, y; intensity_threshold=intensity_thresh)
                end
                peak_results[i] = (pk_mz, pk_int)
                on_stage(:peaks; idx=spec_idx, mz=pk_mz, intensity=pk_int)
            end
            processed_spectra = peak_results # Now contains peak lists
        
        elseif step === :align
            # Alignment requires a reference peak list, typically from the first spectrum
            reference_peaks = isempty(processed_spectra) ? nothing : processed_spectra[1][1]
            if reference_peaks === nothing
                @error "Alignment requires a :peaks step first."; continue
            end
            tol = get(params, :align_tolerance, 0.002)
            for i in 2:length(processed_spectra)
                tgt_peaks, intens = processed_spectra[i]
                warp_func = align_peaks_lowess(reference_peaks, tgt_peaks; tolerance=tol)
                processed_spectra[i] = (warp_func(tgt_peaks), intens)
                on_stage(:align; idx=indices[i], mz=processed_spectra[i][1], intensity=processed_spectra[i][2])
            end

        elseif step === :bin
            all_pks  = [s[1] for s in processed_spectra]
            all_ints = [s[2] for s in processed_spectra]
            tol  = get(params, :bin_tolerance, 0.002)
            freq = get(params, :bin_min_frequency, 0.25)
            mat, mz_bins = bin_peaks(all_pks, all_ints, tol; frequency_threshold=freq)
            final_result = FeatureMatrix(mat, mz_bins, indices)
            break # Binning is the last step
        end
    end

    if final_result !== nothing
        return final_result
    else
        @warn "Pipeline finished without a :bin step; returning processed spectra."
        return processed_spectra
    end
end

# =============================================================================
# 9) Quality Control Metrics
# =============================================================================

"""
    calculate_ppm_error(measured_mz::Float64, theoretical_mz::Float64) -> Float64

Calculates mass accuracy in parts-per-million (PPM).

# Formula
PPM = 10⁶ × |measured_mz - theoretical_mz| / theoretical_mz
"""
function calculate_ppm_error(measured_mz::Real, theoretical_mz::Real)
    if theoretical_mz == 0
        return Inf
    end
    return 1e6 * abs(Float64(measured_mz) - Float64(theoretical_mz)) / Float64(theoretical_mz)
end

"""
    calculate_ppm_error_bulk(measured_mz::Vector{Float64}, theoretical_mz::Vector{Float64}) -> Vector{Float64}

Calculates PPM errors for multiple mass values.
"""
function calculate_ppm_error_bulk(measured_mz::Vector{Real}, theoretical_mz::Vector{Real})
    return [calculate_ppm_error(m, t) for (m, t) in zip(measured_mz, theoretical_mz)]
end

"""
    calculate_resolution_fwhm(mz::Float64, profile_mz::Vector{Float64}, 
                            profile_intensity::Vector{Float64}) -> Float64

Calculates mass resolution using Full Width at Half Maximum (FWHM).

# Formula
Resolution = m / Δm, where Δm is FWHM

# Arguments
- `mz`: Peak centroid m/z
- `profile_mz`: Full m/z array from profile data
- `profile_intensity`: Full intensity array from profile data

# Returns
Resolution or NaN if cannot be calculated
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
    analyze_mass_accuracy(msi_data, reference_peaks; ppm_tolerance=20.0)

Analyzes mass accuracy across the dataset using known reference peaks.

# Arguments
- `msi_data`: Your MSI dataset
- `reference_peaks`: Dict of theoretical m/z values -> compound names
- `ppm_tolerance`: Initial tolerance for peak matching

# Returns
Comprehensive mass accuracy report
"""
function analyze_mass_accuracy(msi_data, reference_peaks::Dict{Float64,String}; 
                             ppm_tolerance::Float64=5.0, sample_spectra=100)
    println("\n[ MASS ACCURACY ANALYSIS ]")
    println("PPM tolerance: ", ppm_tolerance, " ppm")
    println("Spectra to sample: ", sample_spectra)
    
    theoretical_mz = sort(collect(keys(reference_peaks)))
    ppm_errors = Float64[]
    matched_peaks = Tuple{Float64,Float64,String}[]  # (theoretical, measured, compound)
    
    # Sample spectra across the dataset
    if sample_spectra >= length(msi_data.spectra_metadata)
        spectrum_indices = 1:length(msi_data.spectra_metadata)
    else
        spectrum_indices = round.(Int, range(1, length(msi_data.spectra_metadata), length=sample_spectra))
    end
    
    for idx in spectrum_indices
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
        return (mean_ppm=NaN, std_ppm=NaN, min_ppm=NaN, max_ppm=NaN, optimal_ppm=NaN, n_matches=0, matched_peaks=[], all_ppm_errors=[])
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
    get_common_calibration_standards(standard_type::Symbol)

Returns common calibration masses for different instrument types.

# Supported standards
- `:maldi_pos`: Common MALDI-TOF positive mode calibrants
- `:maldi_neg`: Common MALDI-TOF negative mode calibrants  
- `:esi_pos`: ESI positive mode calibrants
- `:lcms`: LC-MS commonly used standards
"""
function get_common_calibration_standards(standard_type::Symbol=:maldi_pos)
    standards = Dict{Float64,String}()
    
    if standard_type == :maldi_pos
        standards = Dict(
            104.10754 => "C5H4N2 (Imidazole)",
            175.11995 => "C6H15O4P (Glycerophosphocholine fragment)",
            226.15687 => "C10H20NO4P (Phosphocholine)",
            322.04810 => "[Glu1]-Fibrinopeptide B fragment",
            379.09247 => "C12H22O11 (Sucrose)",
            515.32539 => "C26H52NO7P (PC(16:0/0:0))",
            622.02896 => "C20H12O5S2 (1-Hydroxypyrene-3,6,8-trisulfate)",
            757.39917 => "C37H74NO8P (PC(34:1))",
            1046.54198 => "Angiotensin I",
            1296.68477 => "ACTH clip 1-17",
            1570.67744 => "ACTH clip 18-39",
            2465.19829 => "ACTH clip 7-38"
        )
    elseif standard_type == :maldi_neg
        standards = Dict(
            112.98563 => "C2F3O2 (Trifluoroacetate)",
            152.99568 => "C2F6S (Perfluoroethylsulfonate)",
            214.00166 => "C4F7O2 (Heptafluorobutyrate)",
            264.93278 => "C6F6 (Hexafluorobenzene)",
            362.96198 => "C8F15O2 (Perfluorooctanoate)",
            466.96714 => "C10F17O2S (Perfluorooctanesulfonate)"
        )
    elseif standard_type == :esi_pos
        standards = Dict(
            118.08626 => "C5H12NO2 (Valine)",
            175.11900 => "C6H15O4P (Phosphocholine fragment)", 
            524.26496 => "C23H48NO7P (LysoPC(16:0))",
            622.02896 => "C20H12O5S2 (Standard)",
            922.00980 => "C18H18O6N3S3 (Ultramark 1621)"
        )
    end
    
    return standards
end

"""
    generate_qc_report(msi_data; reference_peaks, output_dir)

Generates a comprehensive QC report including mass accuracy and resolution.
"""
function generate_qc_report(msi_data, filename::String; reference_peaks=nothing, output_dir="qc_results", sample_spectra=100)
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
    accuracy_report = analyze_mass_accuracy(msi_data, reference_peaks, sample_spectra=sample_spectra)
    
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

