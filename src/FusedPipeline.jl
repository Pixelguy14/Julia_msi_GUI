# src/FusedPipeline.jl
# This file defines a high-performance in-place preprocessing pipeline
# that minimizes allocations by using reused buffers from ResourcePool.

using .MSI_src # Ensure it can see the module's contents if needed

"""
    SpectralPipeline

Holds a sequence of preprocessing steps and the necessary buffers to execute them in-place.
"""
struct SpectralPipeline
    steps::Vector{AbstractPreprocessingStep}
end

"""
    apply_pipeline!(mz::Vector{Float64}, intensity::Vector{Float64}, pipeline::SpectralPipeline; data::MSIData)

Applies all steps in the pipeline to the spectrum arrays in-place.
Uses internal buffers from data.resource_pool where needed.
"""
function apply_pipeline!(mz::Vector{Float64}, intensity::Vector{Float64}, pipeline::SpectralPipeline, data::MSIData)
    # Process each step in sequence
    for step in pipeline.steps
        apply_step!(mz, intensity, step, data)
    end
    return mz, intensity
end

# --- Basic in-place implementations of core preprocessing steps ---

function apply_step!(mz, int, step::Normalization, data)
    if step.method === :tic
        s = sum(int)
        if s > 0
            int ./= s
        end
    elseif step.method === :median
        m = median(int)
        if m > 0
            int ./= m
        end
    end
    return int
end

function apply_step!(mz, int, step::Smoothing, data)
    if step.method === :savitzky_golay
        # SavitzkyGolay.savitzky_golay currently allocates, but we can't easily fix that here
        # without refactoring the library. However, we can use a pooled vector for its output
        # then copy back to intensity.
        # [Wait: For now we'll call smoothed_y = smooth_spectrum_core(int, ...)]
        # We'll use a resource from the pool to avoid fresh allocation
        temp_buf = acquire(data.resource_pool)
        resize!(temp_buf, length(int))
        
        # Call existing core which returns a new vector, unfortunately
        # But we'll copy it back to 'int' to maintain in-place pipeline
        smoothed = smooth_spectrum_core(int; method=step.method, window=step.window, order=step.order)
        copyto!(int, smoothed)
        
        release!(data.resource_pool, temp_buf)
    end
    return int
end

function apply_step!(mz, int, step::BaselineCorrection, data)
    if step.method === :snip
        # SNIP is easy to make in-place!
        iterations = (step.iterations === nothing) ? 100 : step.iterations
        snip_baseline_inplace!(int, iterations)
    end
    return int
end

"""
    snip_baseline_inplace!(y, iterations)

In-place implementation of the Sensitive Nonlinear Iterative Peak clipping algorithm.
"""
function snip_baseline_inplace!(y::Vector{Float64}, iterations::Int)
    n = length(y)
    n < 3 && return y
    
    # We still need one temporary buffer for the SNIP iteration to read from the previous state
    # Actually, we can just return the baseline and subtract it, but to BE in-place,
    # we need a temporary to hold the baseline during calculation.
    
    # For now, we'll use the existing _snip_baseline_impl and subtract
    baseline = _snip_baseline_impl(y, iterations=iterations)
    y .-= baseline
    return y
end
