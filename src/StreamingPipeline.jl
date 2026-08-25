# src/StreamingPipeline.jl
# ============================================================================
# The Streaming Pipeline Executor
#
# This module provides `process_dataset!`, the Sprint 2 master function that
# streams spectral data through an in-place kernel chain and accumulates
# results into a SparseMatrixCSC without ever holding more than 1 spectrum
# per thread in RAM.
#
# Architecture:
#   1. _iterate_spectra_fast → Mmap zero-copy views
#   2. copyto!(writable_buf, view) → makes mutable copy for kernels
#   3. Kernel chain: smooth! → baseline! → peaks → bin
#   4. Thread-local (I, J, V) sparse accumulators
#   5. Final sparse(I, J, V, num_bins, num_spectra) assembly
#
# This works alongside the existing execute_full_preprocessing in
# PreprocessingPipeline.jl — it does NOT replace the app.jl integration.
# ============================================================================

using SparseArrays
using Printf

# =============================================================================
# Configuration Structs
# =============================================================================

"""
    StreamingStep

Represents a single step in the streaming pipeline.
"""
struct StreamingStep
    name::Symbol
    params::Dict{Symbol, Any}
end

"""
    PipelineConfig

Holds the complete configuration for a streaming pipeline execution.

# Fields
- `steps::Vector{StreamingStep}` — ordered sequence of processing steps
- `reference_peaks::Vector{Float64}` — fixed m/z values for calibration (Category B)
- `num_bins::Int` — number of bins for the output feature matrix
- `min_peaks_per_bin::Int` — minimum peak count to keep a bin
- `frequency_threshold::Float64` — minimum fraction of spectra a bin must appear in (0.0-1.0)

# Example
```julia
config = PipelineConfig(
    steps = [
        StreamingStep(:smoothing, Dict(:method => :savitzky_golay, :window => 9, :order => 2)),
        StreamingStep(:baseline_correction, Dict(:method => :snip, :iterations => 100)),
        StreamingStep(:normalization, Dict(:method => :tic)),
        StreamingStep(:peak_picking, Dict(:method => :profile, :snr_threshold => 3.0)),
    ],
    num_bins = 2000
)
```
"""
struct PipelineConfig
    steps::Vector{StreamingStep}
    reference_peaks::Vector{Float64}
    num_bins::Int
    min_peaks_per_bin::Int
    frequency_threshold::Float64
end

# Convenience constructor with defaults
function PipelineConfig(; steps::Vector{StreamingStep}=StreamingStep[],
                         reference_peaks::Vector{Float64}=Float64[],
                         num_bins::Int=2000,
                         min_peaks_per_bin::Int=3,
                         frequency_threshold::Float64=0.0)
    return PipelineConfig(steps, reference_peaks, num_bins, min_peaks_per_bin, frequency_threshold)
end

# =============================================================================
# Sparse Accumulator (Thread-Local)
# =============================================================================

"""
    SparseAccumulator

Thread-local accumulator for sparse matrix construction.
Collects (row, col, val) triplets that will be assembled into
a SparseMatrixCSC at the end of the pipeline.
"""
mutable struct SparseAccumulator
    I::Vector{Int}    # Row indices (bin indices)
    J::Vector{Int}    # Column indices (spectrum indices)
    V::Vector{Float64} # Values (intensities)
    lck::Base.Threads.SpinLock

    function SparseAccumulator(capacity_hint::Int=10000)
        acc = new(
            Vector{Int}(undef, 0),
            Vector{Int}(undef, 0),
            Vector{Float64}(undef, 0),
            Base.Threads.SpinLock()
        )
        sizehint!(acc.I, capacity_hint)
        sizehint!(acc.J, capacity_hint)
        sizehint!(acc.V, capacity_hint)
        return acc
    end
end

"""
    accumulate!(acc::SparseAccumulator, spectrum_idx::Int, bin_indices::AbstractVector{Int},
                intensities::AbstractVector{Float64})

Appends peak data for one spectrum into the sparse accumulator.
"""
@inline function accumulate!(acc::SparseAccumulator, spectrum_idx::Int,
                             bin_indices::AbstractVector{Int},
                             intensities::AbstractVector{Float64})
    n = length(bin_indices)
    for k in 1:n
        @inbounds begin
            push!(acc.I, bin_indices[k])
            push!(acc.J, spectrum_idx)
            push!(acc.V, intensities[k])
        end
    end
end

# =============================================================================
# The Pipeline Executor
# =============================================================================

"""
    process_dataset!(data::MSIData, config::PipelineConfig;
                     progress_callback::Union{Function, Nothing}=nothing,
                     masked_indices::Union{AbstractVector{Int}, Nothing}=nothing)

The Sprint 2 master streaming function. Processes an entire MSI dataset through
a kernel chain without holding more than 1 spectrum per thread in RAM.

# Returns
- `SparseMatrixCSC{Float64, Int}`: The feature matrix (bins × spectra)
- `Vector{Float64}`: The m/z bin centers

# Architecture
1. Ensures analytics are computed (for global m/z range)
2. Creates thread-local SparseAccumulators
3. Streams spectra via `_iterate_spectra_fast`
4. Per spectrum: copy view → kernel chain → peak detect → bin → accumulate
5. Merges accumulators → `sparse(I, J, V)`
"""
function process_dataset!(data::MSIData, config::PipelineConfig;
                          progress_callback::Union{Function, Nothing}=nothing,
                          masked_indices::Union{AbstractVector{Int}, Nothing}=nothing)
    
    # --- Step 1: Ensure analytics are computed (provides global m/z range) ---
    if !is_set(data.analytics_ready)
        println("Pre-computing analytics for streaming pipeline...")
        precompute_analytics(data)
    end
    
    # Determine global m/z range for binning
    global_min_mz = Base.Threads.atomic_add!(data.global_min_mz, 0.0)
    global_max_mz = Base.Threads.atomic_add!(data.global_max_mz, 0.0)
    
    if !isfinite(global_min_mz) || !isfinite(global_max_mz) || global_min_mz >= global_max_mz
        @warn "Invalid global m/z range: [$global_min_mz, $global_max_mz]. Cannot bin peaks."
        return spzeros(0, 0), Float64[]
    end
    
    num_bins = config.num_bins
    bin_edges = range(global_min_mz, stop=global_max_mz, length=num_bins + 1)
    bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in 1:num_bins]
    inv_bin_width = 1.0 / step(bin_edges)
    
    num_spectra = length(data.spectra_metadata)
    indices_to_process = masked_indices === nothing ? nothing : masked_indices
    
    # --- Step 2: Create thread-local accumulators ---
    n_threads = Base.Threads.nthreads()
    accumulators = [SparseAccumulator(num_spectra * 10) for _ in 1:n_threads]
    spectra_processed = Base.Threads.Atomic{Int}(0)
    
    # NEW: Create dedicated workspace buffers for each thread.
    # This completely eliminates the need for acquire/release and prevents deadlocks.
    workspaces_mz = [Vector{Float64}(undef, 0) for _ in 1:n_threads]
    workspaces_int = [Vector{Float64}(undef, 0) for _ in 1:n_threads]
    workspaces_scratch = [Vector{Float64}(undef, 0) for _ in 1:n_threads]
    
    # Pre-parse step configuration for fast dispatch in the hot loop
    has_smoothing = false
    has_baseline = false
    has_normalization = false
    has_transform = false
    has_peak_picking = false
    has_calibration = false
    
    smooth_params = Dict{Symbol, Any}()
    baseline_params = Dict{Symbol, Any}()
    norm_params = Dict{Symbol, Any}()
    transform_params = Dict{Symbol, Any}()
    peak_params = Dict{Symbol, Any}()
    
    for s in config.steps
        if s.name === :smoothing
            has_smoothing = true
            smooth_params = s.params
        elseif s.name === :baseline_correction
            has_baseline = true
            baseline_params = s.params
        elseif s.name === :normalization
            has_normalization = true
            norm_params = s.params
        elseif s.name === :stabilization || s.name === :intensity_transformation
            has_transform = true
            transform_params = s.params
        elseif s.name === :peak_picking
            has_peak_picking = true
            peak_params = s.params
        elseif s.name === :calibration
            has_calibration = true
        end
    end
    
    reference_masses = config.reference_peaks
    
    # --- Step 3: Stream and process ---
    start_time = time_ns()
    
    # Use let block to capture all variables cleanly for the closure
    let data=data, accumulators=accumulators, spectra_processed=spectra_processed,
        bin_edges=bin_edges, num_bins=num_bins, inv_bin_width=inv_bin_width,
        global_min_mz=global_min_mz,
        workspaces_mz=workspaces_mz, workspaces_int=workspaces_int, workspaces_scratch=workspaces_scratch,
        has_smoothing=has_smoothing, has_baseline=has_baseline,
        has_normalization=has_normalization, has_transform=has_transform,
        has_peak_picking=has_peak_picking, has_calibration=has_calibration,
        smooth_params=smooth_params, baseline_params=baseline_params,
        norm_params=norm_params, transform_params=transform_params,
        peak_params=peak_params, reference_masses=reference_masses
        
        _iterate_spectra_fast(data, indices_to_process) do idx, mz_view, int_view
            thread_id = Base.Threads.threadid()
            acc = accumulators[thread_id]
            
            # --- Grab Thread-Local Workspaces ---
            # No locking, no blocking, guaranteed to be available
            mz_buf = workspaces_mz[thread_id]
            int_buf = workspaces_int[thread_id]
            scratch_buf = workspaces_scratch[thread_id]
            
            resize!(mz_buf, length(mz_view))
            resize!(int_buf, length(int_view))
            resize!(scratch_buf, length(int_view))
            copyto!(mz_buf, mz_view)
            copyto!(int_buf, int_view)
            
            # --- Kernel Chain (in pipeline order) ---
                
                # Category B: Fixed-reference calibration
                if has_calibration && !isempty(reference_masses)
                    calibrate_inplace!(mz_buf, int_buf, reference_masses)
                end
                
                # Category A: Intensity transformation
                if has_transform
                    transform_inplace!(int_buf, get(transform_params, :method, :sqrt))
                end
                
                # Category A: Smoothing
                if has_smoothing
                    smooth_inplace!(int_buf, scratch_buf, data;
                        method=get(smooth_params, :method, :savitzky_golay),
                        window=get(smooth_params, :window, 9),
                        order=get(smooth_params, :order, 2))
                end
                
                # Category A: Baseline correction
                if has_baseline
                    baseline_subtract_inplace!(int_buf, scratch_buf, data;
                        method=get(baseline_params, :method, :snip),
                        iterations=get(baseline_params, :iterations, 100),
                        window=get(baseline_params, :window, 20))
                end
                
                # Category A: Normalization
                if has_normalization
                    normalize_inplace!(int_buf, get(norm_params, :method, :tic))
                end
                
                # --- Peak Detection & Binning ---
                if has_peak_picking
                    detect_peaks_streaming(mz_buf, int_buf, scratch_buf;
                        method=get(peak_params, :method, :profile),
                        snr_threshold=Float64(get(peak_params, :snr_threshold, 3.0)),
                        half_window=Int(get(peak_params, :half_window, 10)),
                        min_peak_prominence=Float64(get(peak_params, :min_peak_prominence, 0.1)),
                        merge_peaks_tolerance=Float64(get(peak_params, :merge_peaks_tolerance, 0.002))) do peak_mz, peak_int
                        
                        # Bin each discovered peak directly
                        bin_idx = trunc(Int, (peak_mz - global_min_mz) * inv_bin_width) + 1
                        bin_idx = clamp(bin_idx, 1, num_bins)
                        
                        push!(acc.I, bin_idx)
                        push!(acc.J, idx)
                        push!(acc.V, peak_int)
                    end
                else
                    # No peak picking: bin raw intensity directly
                    @inbounds for i in eachindex(mz_buf)
                        bin_idx = trunc(Int, (mz_buf[i] - global_min_mz) * inv_bin_width) + 1
                        bin_idx = clamp(bin_idx, 1, num_bins)
                        
                        push!(acc.I, bin_idx)
                        push!(acc.J, idx)
                        push!(acc.V, int_buf[i])
                    end
                end
                
                Base.Threads.atomic_add!(spectra_processed, 1)
        end
    end
    
    # --- Step 4: Merge thread-local accumulators ---
    total_entries = sum(length(acc.I) for acc in accumulators)
    merged_I = Vector{Int}(undef, total_entries)
    merged_J = Vector{Int}(undef, total_entries)
    merged_V = Vector{Float64}(undef, total_entries)
    
    offset = 0
    for acc in accumulators
        n = length(acc.I)
        if n > 0
            copyto!(merged_I, offset + 1, acc.I, 1, n)
            copyto!(merged_J, offset + 1, acc.J, 1, n)
            copyto!(merged_V, offset + 1, acc.V, 1, n)
            offset += n
        end
    end
    
    # --- Step 5: Assemble sparse matrix ---
    # Use max combiner: when multiple peaks map to the same bin for same spectrum,
    # keep the maximum intensity
    feature_matrix = sparse(merged_I, merged_J, merged_V, num_bins, num_spectra, max)
    
    # --- Step 6: Apply frequency threshold if configured ---
    if config.frequency_threshold > 0.0
        # Count how many spectra have a non-zero value in each bin
        bin_presence = vec(sum(feature_matrix .> 0, dims=2))
        min_count = ceil(Int, config.frequency_threshold * num_spectra)
        keep_bins = findall(bin_presence .>= min_count)
        feature_matrix = feature_matrix[keep_bins, :]
        bin_centers = bin_centers[keep_bins]
    end
    
    duration = (time_ns() - start_time) / 1e9
    n_processed = spectra_processed[]
    n_nonzeros = nnz(feature_matrix)
    sparsity = 1.0 - n_nonzeros / (size(feature_matrix, 1) * size(feature_matrix, 2) + 1)
    
    @printf "Streaming pipeline complete: %d spectra processed in %.2f seconds.\n" n_processed duration
    @printf "Feature matrix: %d bins × %d spectra, %d non-zeros (%.1f%% sparse)\n" size(feature_matrix, 1) size(feature_matrix, 2) n_nonzeros sparsity * 100
    @printf "RAM: %.1f MB (vs %.1f MB dense)\n" (n_nonzeros * 16) / 1e6 (size(feature_matrix, 1) * size(feature_matrix, 2) * 8) / 1e6
    
    if progress_callback !== nothing
        progress_callback(1.0)
    end
    
    return feature_matrix, collect(Float64, bin_centers)
end

"""
    save_sparse_matrix(matrix::SparseMatrixCSC, output_path::String)

Exports a highly optimized SparseMatrixCSC array to disk using the standard
Matrix Market Coordinate format (`.mtx`), guaranteeing no bottleneck or OOM crashes 
for extremely large MS dataset persistence.
"""
function save_sparse_matrix(matrix::SparseMatrixCSC{Float64, Int}, output_path::String)
    m, n = size(matrix)
    nnz_val = nnz(matrix)
    # Use streaming I/O with a large buffer for ultra-fast persistence
    open(output_path, "w") do io
        # Write Matrix Market Header
        write(io, "%%MatrixMarket matrix coordinate real general\n")
        write(io, "$m $n $nnz_val\n")
        
        # Directly extract CSC properties (O(1) memory, zero allocation)
        row_indices = rowvals(matrix)
        values_array = nonzeros(matrix)
        
        @inbounds for filter_j in 1:n
            # nzrange returns the index bounds for non-zero elements in column 'j'
            for idx in nzrange(matrix, filter_j)
                i = row_indices[idx]
                v = values_array[idx]
                write(io, "$i $filter_j $v\n")
            end
        end
    end
end
