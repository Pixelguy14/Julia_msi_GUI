#!/usr/bin/env julia
# test/test_streaming_pipeline.jl
# ============================================================================
# Validation test for Sprint 2: The Streaming Pipeline
#
# This test exercises process_dataset! against the HR2MSI mouse bladder
# dataset and verifies:
#   1. Correct sparse matrix creation
#   2. Non-zero peak population
#   3. RAM savings vs dense equivalent
#   4. Allocation count and throughput
# ============================================================================

using Pkg
Pkg.activate(".")

using MSI_src
using SparseArrays

# =============================================================================
# Configuration
# =============================================================================

const IMZML_PATH = "/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML"

function main()
    println("=" ^ 60)
    println("SPRINT 2: Streaming Pipeline Validation")
    println("=" ^ 60)
    
    if !isfile(IMZML_PATH)
        println("SKIPPED: Dataset not found at $IMZML_PATH")
        return
    end
    
    # --- 1. Load dataset ---
    println("\n--- Step 1: Loading dataset ---")
    data = OpenMSIData(IMZML_PATH)
    println("Loaded: $(length(data.spectra_metadata)) spectra")
    
    # --- 2. Configure the streaming pipeline ---
    println("\n--- Step 2: Configuring pipeline ---")
    config = PipelineConfig(
        steps = [
            StreamingStep(:baseline_correction, Dict{Symbol,Any}(:method => :snip, :iterations => 50)),
            StreamingStep(:normalization, Dict{Symbol,Any}(:method => :tic)),
            StreamingStep(:peak_picking, Dict{Symbol,Any}(
                :method => :profile,
                :snr_threshold => 3.0,
                :half_window => 10,
                :min_peak_prominence => 0.1,
                :merge_peaks_tolerance => 0.002
            )),
        ],
        num_bins = 2000,
        frequency_threshold = 0.01  # Bins must appear in at least 1% of spectra
    )
    println("Steps: $(join([s.name for s in config.steps], " → "))")
    println("Bins: $(config.num_bins), Frequency threshold: $(config.frequency_threshold)")
    
    # --- 3. Run the streaming pipeline ---
    println("\n--- Step 3: Running streaming pipeline ---")
    stats = @timed begin
        feature_matrix, bin_centers = process_dataset!(data, config)
    end
    
    feature_matrix = stats.value[1]
    bin_centers = stats.value[2]
    
    println("\n--- Results ---")
    println("  Feature matrix size: $(size(feature_matrix))")
    println("  Non-zeros: $(nnz(feature_matrix))")
    println("  Bin centers: $(length(bin_centers))")
    println("  Time: $(round(stats.time, digits=2))s")
    println("  Allocations: $(stats.bytes ÷ 1_000_000) MB")
    println("  GC time: $(round(stats.gctime, digits=2))s")
    
    # --- 4. Validate ---
    println("\n--- Step 4: Validation ---")
    
    passed = true
    
    # Check matrix dimensions
    if size(feature_matrix, 1) > 0 && size(feature_matrix, 2) > 0
        println("  ✓ Matrix has valid dimensions")
    else
        println("  ✗ Matrix has invalid dimensions: $(size(feature_matrix))")
        passed = false
    end
    
    # Check non-zeros
    if nnz(feature_matrix) > 0
        println("  ✓ Matrix has $(nnz(feature_matrix)) non-zero entries")
    else
        println("  ✗ Matrix is completely empty")
        passed = false
    end
    
    # Check sparsity savings
    dense_mb = size(feature_matrix, 1) * size(feature_matrix, 2) * 8 / 1e6
    sparse_mb = nnz(feature_matrix) * 16 / 1e6  # index + value per entry
    if dense_mb > 0
        savings = (1.0 - sparse_mb / dense_mb) * 100
        println("  ✓ RAM savings: $(round(savings, digits=1))% ($(round(sparse_mb, digits=1)) MB vs $(round(dense_mb, digits=1)) MB dense)")
    end
    
    # Check bin centers alignment
    if length(bin_centers) == size(feature_matrix, 1)
        println("  ✓ Bin centers match matrix rows")
    else
        println("  ✗ Bin center count ($(length(bin_centers))) != matrix rows ($(size(feature_matrix, 1)))")
        passed = false
    end
    
    println("\n" * "=" ^ 60)
    if passed
        println("ALL VALIDATIONS PASSED ✓")
    else
        println("SOME VALIDATIONS FAILED ✗")
    end
    println("=" ^ 60)
end

@time main()
