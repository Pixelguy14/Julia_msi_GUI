using BenchmarkTools
using MSI_src
using Statistics
using DataFrames

# ===================================================================
# HIGH-PRECISION COMPARATIVE SUITE
# ===================================================================

function run_advanced_benchmark(path, mz, tol)
    println("\n" * "="^40)
    println("TARGET: $(basename(path))")
    println("="^40)

    # 1. NEW LIBRARY: Metadata Load (The "Control Tower" startup)
    # This measures how fast the Mmap and Cache system works
    t_load_new = @belapsed OpenMSIData($path)
    
    # 2. NEW LIBRARY: Slice Generation (The "Streaming" speed)
    msi_new = OpenMSIData(path)
    # We use @benchmark to get a distribution (min, mean, max)
    b_slice_new = @benchmark get_mz_slice($msi_new, $mz, $tol)

    # --- Metrics Table ---
    results = DataFrame(
        Metric = ["Metadata Load", "Slice Gen (Min)", "Slice Gen (Mean)", "Allocations"],
        JuliaMSI = [
            "$(round(t_load_new * 1000, digits=2)) ms",
            "$(round(minimum(b_slice_new.times)/1e6, digits=2)) ms",
            "$(round(mean(b_slice_new.times)/1e6, digits=2)) ms",
            "$(b_slice_new.allocs) allocs"
        ]
    )
    
    println(results)
    
    # --- The "Throughput" Test ---
    # How many slices per second can we handle?
    throughput_new = 1.0 / mean(b_slice_new.times/1e9)
    println("\nThroughput: $(round(throughput_new, digits=1)) slices/sec")
    
    return results
end

# Example Run
@time run_advanced_benchmark("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML", 716.053, 0.1)

# For multithread: julia --threads auto --project=. test/new_benchmark_mmap.jl