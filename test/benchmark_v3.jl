# test/benchmark_v3.jl

# ===================================================================
# High-Precision Performance Benchmark Suite (v3) - Cross-Engine
# ===================================================================
# Compares JuliaMSI (New) vs julia_mzML_imzML (Legacy)
# Integrates formal academic metrics, error bars, and log-scale plots.
# ===================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using DataFrames
using CSV
using CairoMakie
using Statistics
using MSI_src

struct BenchmarkCase
    filepath::String
    mz_value::Float64
    mz_tolerance::Float64
    name::String
end

const RESULTS_DIR = joinpath(@__DIR__, "results")
const NUM_REPETITIONS = 10

const BENCHMARK_CASES = [
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_LA-ESI/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML", 436.04, 0.1, "Thaliana Leaf"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Chilli/ltpmsi-chilli.imzML", 420.0, 0.1, "Chilli Pepper"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_DESI/40TopL,10TopR,30BottomL,20BottomR/40TopL,10TopR,30BottomL,20BottomR-centroid.imzML", 885.5, 0.1, "Colon Cancer Human"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML", 716.053, 0.1, "Mouse Urinary Bladder"),
    #BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Leafs/CE1_Leaf_R3.imzML", 306.1, 0.1, "Leaf"),
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
        
        load_stats_new = @timed OpenMSIData(case.filepath)
        msi_data = load_stats_new.value
        load_time_new_s = load_stats_new.time
        mem_load_new_mb = load_stats_new.bytes / 1024^2
        
        # Benchmark slice generation and capture variability
        b_slice_new = @benchmark get_mz_slice($msi_data, $(case.mz_value), $(case.mz_tolerance)) samples=NUM_REPETITIONS seconds=10
        times_new_s = b_slice_new.times ./ 1e9
        mean_time_new_s = mean(times_new_s)
        std_time_new_s = std(times_new_s)
        
        total_time_10_new = load_time_new_s + (10 * mean_time_new_s)
        amortized_sps_new = 10.0 / total_time_10_new
        
        close(msi_data)
        flush_memory()

        # ------------------------------------------------------------
        # 2. julia_mzML_imzML (Legacy Architecture)
        # ------------------------------------------------------------
        println("[julia_mzML_imzML - Legacy]")
        
        mem_load_old_mb = NaN
        load_time_old_s = Inf
        mean_time_old_s = Inf
        std_time_old_s = 0.0
        amortized_sps_old = 0.0
        legacy_status = "Success"
        
        # Skip files > 3 GB in legacy to prevent OOM / SIGKILL
        if file_size > 5000
            println("Skipping legacy benchmark for $(case.name) to prevent SIGKILL OOM.")
            legacy_status = "OOM / Skipped"
        else
            temp_legacy_csv = joinpath(RESULTS_DIR, "temp_legacy_v3.csv")
            rm(temp_legacy_csv, force=true)
            
            legacy_cmd_str = """
            using julia_mzML_imzML, Statistics, BenchmarkTools, CSV, DataFrames
            try
                load_stats = @timed LoadImzml("$(case.filepath)")
                old_data = load_stats.value
                mem_load_mb = load_stats.bytes / 1024^2
                load_time_s = load_stats.time
                
                b_slice = @benchmark GetSlice(\$old_data, $(case.mz_value), $(case.mz_tolerance)) samples=$NUM_REPETITIONS seconds=10
                times_s = b_slice.times ./ 1e9
                
                df = DataFrame(
                    LoadMem_MB=[mem_load_mb], 
                    LoadTime_s=[load_time_s], 
                    SliceMean_s=[mean(times_s)], 
                    SliceStd_s=[std(times_s)]
                )
                CSV.write("$(temp_legacy_csv)", df)
            catch e
                df = DataFrame(LoadMem_MB=[NaN], LoadTime_s=[Inf], SliceMean_s=[NaN], SliceStd_s=[0.0])
                CSV.write("$(temp_legacy_csv)", df)
            end
            """
            
            try
                run(`julia --project=. -e $(legacy_cmd_str)`)
                if isfile(temp_legacy_csv)
                    legacy_df = CSV.read(temp_legacy_csv, DataFrame)
                    mem_load_old_mb = legacy_df.LoadMem_MB[1]
                    load_time_old_s = legacy_df.LoadTime_s[1]
                    mean_time_old_s = legacy_df.SliceMean_s[1]
                    std_time_old_s = legacy_df.SliceStd_s[1]
                    
                    if isfinite(mean_time_old_s)
                        total_time_10_old = load_time_old_s + (10 * mean_time_old_s)
                        amortized_sps_old = 10.0 / total_time_10_old
                    else
                        legacy_status = "Failed"
                    end
                    rm(temp_legacy_csv, force=true)
                else
                    legacy_status = "OOM / Failed"
                end
            catch e
                @warn "Legacy process crashed: $e"
                legacy_status = "OOM / Failed"
            end
        end
        
        flush_memory()

        # ------------------------------------------------------------
        # Derived Metrics Calculation
        # ------------------------------------------------------------
        # Throughput (Slices per second)
        tput_new = 1.0 / mean_time_new_s
        tput_old = isfinite(mean_time_old_s) ? (1.0 / mean_time_old_s) : 0.0
        
        # Memory Efficiency (MB RAM / GB File)
        gb_file = file_size / 1024.0
        mem_eff_new = mem_load_new_mb / gb_file
        mem_eff_old = isfinite(mem_load_old_mb) ? (mem_load_old_mb / gb_file) : NaN
        
        # Fold Speedup (Legacy / JuliaMSI)
        speedup = isfinite(mean_time_old_s) ? (mean_time_old_s / mean_time_new_s) : NaN
        
        # Load Bandwidth (MB/s)
        bandwidth_new = file_size / load_time_new_s
        bandwidth_old = isfinite(load_time_old_s) ? (file_size / load_time_old_s) : 0.0

        push!(results, (
            Dataset = case.name,
            FileSize_MB = file_size,
            
            # New Engine (JuliaMSI)
            New_LoadTime_s = load_time_new_s,
            New_LoadMem_MB = mem_load_new_mb,
            New_SliceMean_s = mean_time_new_s,
            New_SliceStd_s = std_time_new_s,
            New_Throughput_SPS = tput_new,
            New_MemEff_MB_GB = mem_eff_new,
            New_Bandwidth_MBs = bandwidth_new,
            New_Amortized_SPS_10 = amortized_sps_new,
            
            # Legacy Engine
            Legacy_LoadTime_s = load_time_old_s,
            Legacy_LoadMem_MB = mem_load_old_mb,
            Legacy_SliceMean_s = mean_time_old_s,
            Legacy_SliceStd_s = std_time_old_s,
            Legacy_Throughput_SPS = tput_old,
            Legacy_MemEff_MB_GB = mem_eff_old,
            Legacy_Bandwidth_MBs = bandwidth_old,
            Legacy_Amortized_SPS_10 = amortized_sps_old,
            
            # Comparative Deltas
            Fold_Speedup = speedup,
            RAM_Reduction_Pct = isfinite(mem_load_old_mb) ? ((mem_load_old_mb - mem_load_new_mb) / mem_load_old_mb)*100 : NaN,
            Legacy_Status = legacy_status
        ))
    end

    csv_path = joinpath(RESULTS_DIR, "v3_enterprise_benchmarks.csv")
    CSV.write(csv_path, results)
    println("\nData saved to $csv_path")
    
    plot_v3_results(results)
    return results
end

function plot_v3_results(df::DataFrame)
    if isempty(df) return end

    fig = Figure(size=(1800, 1200), fontsize=22)
    x_pos = 1:nrow(df)
    labels = df.Dataset

    # ---------- Panel A: Initial Load Time (Log Scale) ----------
    ax1 = Axis(fig[1, 1], 
               title="Panel A: Metadata/File Load Time", 
               ylabel="Time (seconds) - Log Scale",
               yscale=log10,
               xticks=(x_pos, labels), xticklabelrotation=π/8)
    
    barplot!(ax1, x_pos .- 0.2, df.New_LoadTime_s, color="#10b981", width=0.4, label="JuliaMSI")
    
    # Only plot legacy if it succeeded
    legacy_times = [row.Legacy_Status == "Success" ? row.Legacy_LoadTime_s : 1e-3 for row in eachrow(df)]
    barplot!(ax1, x_pos .+ 0.2, legacy_times, color="#ef4444", width=0.4, label="Legacy")
    axislegend(ax1, position=:lt)

    # ---------- Panel B: Initial RAM Cost (Log Scale) ----------
    ax2 = Axis(fig[1, 2], 
               title="Panel B: Memory Footprint (RAM)", 
               ylabel="Memory (MB) - Log Scale",
               yscale=log10,
               xticks=(x_pos, labels), xticklabelrotation=π/8)
    
    barplot!(ax2, x_pos .- 0.2, df.New_LoadMem_MB, color="#10b981", width=0.4, label="JuliaMSI (Mmap)")
    
    legacy_mems = [row.Legacy_Status == "Success" ? row.Legacy_LoadMem_MB : 1e-1 for row in eachrow(df)]
    barplot!(ax2, x_pos .+ 0.2, legacy_mems, color="#ef4444", width=0.4, label="Legacy (Dense)")
    axislegend(ax2, position=:lt)

    # ---------- Panel C: Slice Throughput & Error Bars ----------
    ax3 = Axis(fig[2, 1:2], 
               title="Panel C: Slice Generation Time (Mean ± SD)", 
               ylabel="Time (seconds) - Log Scale",
               yscale=log10,
               xticks=(x_pos, labels), xticklabelrotation=π/8)
               
    # Plot JuliaMSI times and error bars
    barplot!(ax3, x_pos .- 0.2, df.New_SliceMean_s, color="#10b981", width=0.4, label="JuliaMSI")
    errorbars!(ax3, x_pos .- 0.2, df.New_SliceMean_s, df.New_SliceStd_s, color=:black, whiskerwidth=15)
    
    # Plot Legacy times and error bars (where successful)
    legacy_means = [row.Legacy_Status == "Success" ? row.Legacy_SliceMean_s : 1e-4 for row in eachrow(df)]
    legacy_stds = [row.Legacy_Status == "Success" ? row.Legacy_SliceStd_s : 0.0 for row in eachrow(df)]
    barplot!(ax3, x_pos .+ 0.2, legacy_means, color="#ef4444", width=0.4, label="Legacy")
    errorbars!(ax3, x_pos .+ 0.2, legacy_means, legacy_stds, color=:black, whiskerwidth=15)
    axislegend(ax3, position=:rt)

    save(joinpath(RESULTS_DIR, "v3_enterprise_dashboard.png"), fig)
    println("\nDashboards saved successfully to $(joinpath(RESULTS_DIR, "v3_enterprise_dashboard.png"))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_v3_benchmarks()
end
