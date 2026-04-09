# src/StreamingKernels.jl
# ============================================================================
# In-Place Spectral Kernels for the Streaming Pipeline
# 
# These functions operate on raw (mz, intensity) views from the Sprint 1
# Mmap engine. They write results back to the input buffers using .= to
# achieve zero-allocation processing per spectrum.
#
# Design contract:
#   - All !-suffixed functions modify their arguments in-place
#   - If a kernel needs temporary storage, it borrows from data.resource_pool
#   - No function creates MutableSpectrum objects
# ============================================================================

using Statistics: mean, median

# =============================================================================
# Category A: Purely Streamable Kernels
# =============================================================================

"""
    normalize_inplace!(intensity::AbstractVector{<:Real}, method::Symbol)

Normalizes intensity values in-place. Supports :tic, :median, :rms.
Zero-allocation for the normalization itself.
"""
@inline function normalize_inplace!(intensity::AbstractVector{<:Real}, method::Symbol)
    if method === :tic
        s = sum(intensity)
        if s > 0
            intensity ./= s
        end
    elseif method === :median
        m = median(intensity)
        if m > 0
            intensity ./= m
        end
    elseif method === :rms
        s = sqrt(sum(abs2, intensity) / length(intensity))
        if s > 0
            intensity ./= s
        end
    end
    return intensity
end

"""
    transform_inplace!(intensity::AbstractVector{Float64}, method::Symbol)

Applies intensity transformation in-place. Supports :sqrt, :log1p, :log, :log2, :log10.
"""
@inline function transform_inplace!(intensity::AbstractVector{Float64}, method::Symbol)
    if method === :sqrt
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = sqrt(max(0.0, intensity[i]))
        end
    elseif method === :log1p
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = log1p(max(0.0, intensity[i]))
        end
    elseif method === :log
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = log(max(eps(Float64), intensity[i]))
        end
    elseif method === :log2
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = log2(max(eps(Float64), intensity[i]))
        end
    elseif method === :log10
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = log10(max(eps(Float64), intensity[i]))
        end
    end
    return intensity
end

"""
    smooth_inplace!(intensity::AbstractVector{Float64}, data::MSIData;
                    method::Symbol=:savitzky_golay, window::Int=9, order::Int=2)

Smooths intensity in-place using a temporary buffer from the resource pool.
The SavitzkyGolay library allocates internally, but we copy the result back
to the original buffer and return the pool buffer.
"""
function smooth_inplace!(intensity::AbstractVector{Float64}, scratch::AbstractVector{Float64}, data::MSIData;
                         method::Symbol=:savitzky_golay, window::Int=9, order::Int=2)
    n = length(intensity)
    if n < 3
        return intensity
    end

    if method === :savitzky_golay
        win = isodd(window) ? window : window + 1
        if n < win
            return intensity
        end
        # SavitzkyGolay handles its own math but causes mild allocation.
        res = SavitzkyGolay.savitzky_golay(collect(intensity), win, order)
        @inbounds for i in eachindex(intensity)
            intensity[i] = max(0.0, res.y[i])
        end
    elseif method === :moving_average
        copyto!(scratch, intensity)
        
        half_w = div(window, 2)
        @inbounds for i in 1:n
            s_idx = max(1, i - half_w)
            e_idx = min(n, i + half_w)
            s = 0.0
            @simd for j in s_idx:e_idx
                s += scratch[j]
            end
            intensity[i] = max(0.0, s / (e_idx - s_idx + 1))
        end
    end
    return intensity
end

"""
    baseline_subtract_inplace!(intensity::AbstractVector{Float64}, data::MSIData;
                               method::Symbol=:snip, iterations::Int=100, window::Int=20)

Subtracts baseline from intensity in-place. Uses two pool buffers for the
SNIP ping-pong iteration to avoid any heap allocation in the hot loop.
"""
function baseline_subtract_inplace!(intensity::AbstractVector{Float64}, scratch::AbstractVector{Float64}, data::MSIData;
                                    method::Symbol=:snip, iterations::Int=100, window::Int=20)
    n = length(intensity)
    if n < 3
        return intensity
    end
    
    if method === :snip
        copyto!(scratch, intensity)
        
        for k in 1:iterations
            prev_val = scratch[1]
            scratch[1] = min(scratch[1], scratch[2])
            
            @inbounds for i in 2:n-1
                curr_val = scratch[i]
                scratch[i] = min(curr_val, 0.5 * (prev_val + scratch[i+1]))
                prev_val = curr_val
            end
            scratch[n] = min(scratch[n], prev_val)
        end
        
        @inbounds @simd for i in 1:n
            intensity[i] = max(0.0, intensity[i] - scratch[i])
        end
    elseif method === :convex_hull
        baseline = convex_hull_baseline(intensity)
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = max(0.0, intensity[i] - baseline[i])
        end
    elseif method === :median
        baseline = median_baseline(intensity; window=window)
        @inbounds @simd for i in eachindex(intensity)
            intensity[i] = max(0.0, intensity[i] - baseline[i])
        end
    end
    return intensity
end

"""
    detect_peaks_streaming(mz::AbstractVector, intensity::AbstractVector;
                          method::Symbol=:profile, snr_threshold::Float64=3.0,
                          half_window::Int=10, min_peak_prominence::Float64=0.1,
                          merge_peaks_tolerance::Float64=0.002)

Detects peaks and returns a vector of (mz, intensity) tuples.
This delegates to existing _core functions but returns a lightweight format
suitable for sparse accumulation (no NamedTuple overhead in the hot path).
"""
function detect_peaks_streaming(callback::Function, mz::AbstractVector{Float64}, intensity::AbstractVector{Float64}, scratch::AbstractVector{Float64};
                               method::Symbol=:profile, snr_threshold::Float64=3.0,
                               half_window::Int=10, min_peak_prominence::Float64=0.1,
                               merge_peaks_tolerance::Float64=0.002)
    n = length(intensity)
    if n < 3
        return
    end

    if method === :profile || method === :wavelet
        # Zero-allocation noisy estimation (using mean of bottom half)
        sum_i = 0.0
        @simd for i in 1:n
            sum_i += intensity[i]
        end
        mean_i = sum_i / n
        
        sum_noise = 0.0
        count_noise = 0
        @inbounds for i in 1:n
            if intensity[i] < mean_i
                sum_noise += intensity[i]
                count_noise += 1
            end
        end
        # Use * 1.5 as an approximation to MAD
        noise_level = count_noise > 0 ? (sum_noise / count_noise) * 1.5 + eps(Float64) : mean_i + eps(Float64)
        
        # We will use the scratch buffer to store candidate indices to avoid allocating `Int[]`
        # Because scratch is Float64, we can safely store integer indices up to 2^53 exactly.
        num_candidates = 0
        
        @inbounds for i in 2:n-1
            if intensity[i] > snr_threshold * noise_level
                left = max(1, i - half_window)
                right = min(n, i + half_window)
                
                is_max = true
                for j in left:right
                    if intensity[j] > intensity[i]
                        is_max = false
                        break
                    end
                end
                
                if is_max
                    # Compute prominence
                    min_left = intensity[i]
                    for j in left:i
                        if intensity[j] < min_left
                            min_left = intensity[j]
                        end
                    end
                    min_right = intensity[i]
                    for j in i:right
                        if intensity[j] < min_right
                            min_right = intensity[j]
                        end
                    end
                    
                    prominence = intensity[i] - max(min_left, min_right)
                    
                    if prominence > min_peak_prominence * intensity[i]
                        num_candidates += 1
                        scratch[num_candidates] = i
                    end
                end
            end
        end
        
        # Merge close peaks
        if num_candidates > 0
            if merge_peaks_tolerance > 0
                last_idx = trunc(Int, scratch[1])
                # We emit the first peak lazily down below, so let's compact them in place
                num_merged = 1
                
                for i in 2:num_candidates
                    idx = trunc(Int, scratch[i])
                    if (mz[idx] - mz[last_idx]) > merge_peaks_tolerance
                        num_merged += 1
                        scratch[num_merged] = idx
                        last_idx = idx
                    elseif intensity[idx] > intensity[last_idx]
                        scratch[num_merged] = idx
                        last_idx = idx
                    end
                end
                num_candidates = num_merged
            end
            
            # Emit merged peaks
            for i in 1:num_candidates
                idx = trunc(Int, scratch[i])
                callback(mz[idx], intensity[idx])
            end
        end
        
    elseif method === :centroid
        # Just use any value over snr_threshold * mean_noise
        sum_i = sum(intensity)
        mean_i = sum_i / n
        noise_level = mean_i + eps(Float64)
        
        @inbounds for i in 1:n
            if intensity[i] > snr_threshold * noise_level
                callback(mz[i], intensity[i])
            end
        end
    end
end

# =============================================================================
# Category B: Conditionally Streamable Kernels (Fixed-Reference)
# =============================================================================

"""
    calibrate_inplace!(mz::Vector{Float64}, intensity::AbstractVector,
                       reference_masses::Vector{Float64}; ppm_tolerance::Float64=20.0)

Calibrates the m/z axis in-place using a fixed dictionary of internal standard
reference masses. This is streamable because the reference is constant.

Returns `true` if calibration was applied, `false` if insufficient peaks were found.
"""
function calibrate_inplace!(mz::Vector{Float64}, intensity::AbstractVector,
                            reference_masses::Vector{Float64}; ppm_tolerance::Float64=20.0)
    matched_peaks = find_calibration_peaks_core(mz, intensity, reference_masses; 
                                                 ppm_tolerance=ppm_tolerance)
    if length(matched_peaks) < 2
        return false # Insufficient reference peaks
    end
    
    measured = sort(collect(values(matched_peaks)))
    theoretical = sort(collect(keys(matched_peaks)))
    itp = linear_interpolation(measured, theoretical, extrapolation_bc=Line())
    
    # Apply calibration in-place
    @inbounds for i in eachindex(mz)
        mz[i] = itp(mz[i])
    end
    return true
end
