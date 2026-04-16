# test/benchmark_v3.jl

# ===================================================================
# High-Precision Performance Benchmark Suite (v3)
# ===================================================================
# This script integrates the robust statistical sampling of `BenchmarkTools`
# with the visual comparative mechanics against the legacy library.
# It captures the paradigm shift from "Time per slice" to "Pipeline Throughput".
# ===================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using DataFrames
using CSV
using CairoMakie
using Statistics
using MSI_src
using julia_mzML_imzML

struct BenchmarkCase
    filepath::String
    mz_value::Float64
    mz_tolerance::Float64
    name::String
end

const RESULTS_DIR = joinpath(@__DIR__, "results")

# List to benchmark
const BENCHMARK_CASES = [
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Chilli/ltpmsi-chilli.imzML", 420.0, 0.1, "Chilli Pepper"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_DESI/40TopL,10TopR,30BottomL,20BottomR/40TopL,10TopR,30BottomL,20BottomR-centroid.imzML", 885.5, 0.1, "Colon Cancer Human"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML", 716.053, 0.1, "Mouse Urinary Bladder"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Leafs/CE1_Leaf_R3.imzML", 306.1, 0.1, "Leaf"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Liv2_imzML_TIMSConvert-selected/Liv2.imzML", 796.18, 0.1, "Liver Cut"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/salida_Tims/Stomach_DHB.imzML", 804.3, 0.1, "Mouse Stomach 4GB")
]

function get_total_file_size_mb(filepath::String)
    imzml_size = isfile(filepath) ? filesize(filepath) : 0
    ibd_path = replace(filepath, r"\.(imzML|imzml)$" => ".ibd")
    ibd_size = isfile(ibd_path) ? filesize(ibd_path) : 0
    return round((imzml_size + ibd_size) / 1024^2, digits=2)
end

function flush_memory()
    GC.gc(true)
    if Sys.islinux()
        # Force the OS to reclaim memory from the glibc allocator
        ccall(:malloc_trim, Int32, (Int32,), 0)
    end
end

function run_v3_benchmarks()
    mkpath(RESULTS_DIR)
    results = DataFrame()

    println("="^60)
    println("STARTING ENTERPRISE BENCHMARK SUITE (v3)")
    println("="^60)

    for case in BENCHMARK_CASES
        if !isfile(case.filepath)
            @warn "File not found: $(case.filepath). Skipping."
            continue
        end

        println("\n--- Target: $(case.name) ---")
        file_size = get_total_file_size_mb(case.filepath)

        # ------------------------------------------------------------
        # 1. JuliaMSI (New Architecture)
        # ------------------------------------------------------------
        println("[JuliaMSI - New Engine]")
        flush_memory()
        
        # Load Phase (Metadata Only / Memory Mapping)
        load_stats_new = @timed OpenMSIData(case.filepath)
        msi_data = load_stats_new.value
        load_time_new_s = load_stats_new.time
        mem_load_new_mb = load_stats_new.bytes / 1024^2
        
        # Slicing Phase (High Precision)
        b_slice_new = @benchmark get_mz_slice($msi_data, $(case.mz_value), $(case.mz_tolerance)) samples=10 seconds=5
        mean_time_new_s = mean(b_slice_new.times) / 1e9
        
        # Calculate Amortized Throughput (10 slices)
        # Includes the 'Loading Wall' penalty
        total_time_10_new = load_time_new_s + (10 * mean_time_new_s)
        amortized_sps_new = 10.0 / total_time_10_new
        
        close(msi_data)
        flush_memory()

        # ------------------------------------------------------------
        # 2. julia_mzML_imzML (Legacy Architecture)
        # ------------------------------------------------------------
        println("[julia_mzML_imzML - Legacy]")
        
        load_stats_old = try
            @timed LoadImzml(case.filepath)
        catch e
            @warn "Legacy load failed: $e"
            (time=Inf, bytes=Inf, value=nothing)
        end
        
        old_data = load_stats_old.value
        mem_load_old_mb = load_stats_old.bytes / 1024^2
        
        load_time_old_s = Inf
        mean_time_old_s = Inf
        amortized_sps_old = 0.0
        
        if old_data !== nothing
            b_slice_old = try
                # Legacy can be extremely slow, constrain it heavily
                @benchmark GetSlice($old_data, $(case.mz_value), $(case.mz_tolerance)) samples=3 seconds=10
            catch e
                @warn "Legacy slice failed: $e"
                nothing
            end
            
            if b_slice_old !== nothing
                mean_time_old_s = mean(b_slice_old.times) / 1e9
                
                # Accurately reflect the massive load time in the throughput
                load_time_old_s = load_stats_old.time
                total_time_10_old = load_time_old_s + (10 * mean_time_old_s)
                amortized_sps_old = 10.0 / total_time_10_old
            end
        end
        
        flush_memory()

        # ------------------------------------------------------------
        # Record
        # ------------------------------------------------------------
        push!(results, (
            Dataset = case.name,
            FileSize_MB = file_size,
            
            # Legacy Metrics
            Legacy_LoadMem_MB = mem_load_old_mb,
            Legacy_LoadTime_s = load_stats_old.time,
            Legacy_Amortized_SPS_10 = amortized_sps_old,
            
            # New Metrics
            New_LoadMem_MB = mem_load_new_mb,
            New_LoadTime_s = load_time_new_s,
            New_Amortized_SPS_10 = amortized_sps_new,
            
            # Competitive Deltas
            RAM_Reduction_Pct = isfinite(mem_load_old_mb) ? ((mem_load_old_mb - mem_load_new_mb) / mem_load_old_mb)*100 : NaN,
            UX_Speedup_Factor = isfinite(amortized_sps_old) ? (amortized_sps_new / amortized_sps_old) : NaN
        ))
        
        println("  > RAM Reduction: $(round(results[end, :RAM_Reduction_Pct], digits=2))%")
        println("  > Amortized Throughput (10 Slices): $(round(amortized_sps_old, digits=2)) -> $(round(amortized_sps_new, digits=2)) slices/sec")
    end

    csv_path = joinpath(RESULTS_DIR, "v3_enterprise_benchmarks.csv")
    CSV.write(csv_path, results)
    println("\nData saved to $csv_path")
    
    plot_v3_results(results)
    return results
end

function plot_v3_results(df::DataFrame)
    if isempty(df) return end

    fig = Figure(size=(1800, 1200), fontsize=24)
    x_pos = 1:nrow(df)
    labels = df.Dataset

    # ---------- Plot 1: The RAM Revolution ----------
    ax1 = Axis(fig[1, 1], 
               title="Initial RAM Cost (Loading & Mmap)", 
               ylabel="Memory (MB)",
               xticks=(x_pos, labels), xticklabelrotation=π/8)
    
    barplot!(ax1, x_pos .- 0.2, df.New_LoadMem_MB, color="#10b981", width=0.4, label="JuliaMSI (Mmap Lazy)")
    barplot!(ax1, x_pos .+ 0.2, df.Legacy_LoadMem_MB, color="#ef4444", width=0.4, label="Legacy (Dense Matrix)")
    axislegend(ax1, position=:lt)

    # ---------- Plot 2: Throughput Leap ----------
    ax2 = Axis(fig[1, 2], 
               title="Amortized Throughput (10 Slices - Higher is Better)", 
               ylabel="Slices / Second (Inc. Load Time)",
               xticks=(x_pos, labels), xticklabelrotation=π/8)
               
    barplot!(ax2, x_pos .- 0.2, df.New_Amortized_SPS_10, color="#10b981", width=0.4, label="JuliaMSI")
    barplot!(ax2, x_pos .+ 0.2, df.Legacy_Amortized_SPS_10, color="#ef4444", width=0.4, label="Legacy")
    axislegend(ax2, position=:lt)

    # ---------- Plot 3: The Scaling Wall (File Size vs Throughput) ----------
    ax3 = Axis(fig[2, 1:2], 
               title="Performance Scaling: UX Throughput vs Dataset Size", 
               xlabel="Dataset File Size (MB)", 
               ylabel="Amortized Throughput (Slices/sec)")
               
    scatterlines!(ax3, df.FileSize_MB, df.New_Amortized_SPS_10, color="#10b981", markersize=15, linewidth=4, label="JuliaMSI Engine")
    scatterlines!(ax3, df.FileSize_MB, df.Legacy_Amortized_SPS_10, color="#ef4444", markersize=15, linewidth=4, label="Legacy Engine")
    axislegend(ax3, position=:rt)

    save(joinpath(RESULTS_DIR, "v3_enterprise_dashboard.png"), fig)
    println("\nDashboards saved successfully.")
end

# Ensure we process if run as the main script
if abspath(PROGRAM_FILE) == @__FILE__
    run_v3_benchmarks()
end
