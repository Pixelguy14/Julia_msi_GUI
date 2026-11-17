# test/benchmark_comparison.jl

# ===================================================================
# Cross-Language Performance Benchmark Suite: JuliaMSI vs R (MALDIquant & Cardinal)
# ===================================================================
# Compares JuliaMSI (Julia) vs MALDIquant (R) vs Cardinal (R)
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

function run_cross_benchmarks()
    mkpath(RESULTS_DIR)
    results = DataFrame()

    println("="^60)
    println("STARTING CROSS-LANGUAGE BENCHMARK (JuliaMSI vs R)")
    println("="^60)

    has_r = try
        success(`Rscript --version`)
    catch
        false
    end

    if !has_r
        @warn "Rscript is not available. Please install R first to compare with MALDIquant/Cardinal."
    end

    for case in BENCHMARK_CASES
        if !isfile(case.filepath)
            @warn "File not found: $(case.filepath). Skipping."
            continue
        end

        println("\n--- Target Dataset: $(case.name) ---")
        file_size = get_total_file_size_mb(case.filepath)

        # ------------------------------------------------------------
        # 1. JuliaMSI (New Engine)
        # ------------------------------------------------------------
        println("[JuliaMSI - New Engine]")
        flush_memory()
        
        load_stats_new = @timed OpenMSIData(case.filepath)
        msi_data = load_stats_new.value
        load_time_new_s = load_stats_new.time
        mem_load_new_mb = load_stats_new.bytes / 1024^2
        
        # Benchmark slice generation
        b_slice_new = @benchmark get_mz_slice($msi_data, $(case.mz_value), $(case.mz_tolerance)) samples=NUM_REPETITIONS seconds=10
        times_new_s = b_slice_new.times ./ 1e9
        mean_time_new_s = mean(times_new_s)
        std_time_new_s = std(times_new_s)
        
        close(msi_data)
        flush_memory()

        push!(results, (
            Dataset = case.name,
            FileSize_MB = file_size,
            Library = "JuliaMSI",
            LoadTime_s = load_time_new_s,
            LoadMem_MB = mem_load_new_mb,
            AvgSliceTime_s = mean_time_new_s,
            SliceStd_s = std_time_new_s,
            Status = "Success"
        ))

        # ------------------------------------------------------------
        # 2. Run R MALDIquant (Isolated process)
        # ------------------------------------------------------------
        if has_r
            # Skip MALDIquant for files > 3 GB to prevent system SIGKILL
            if file_size > 5000
                println("Skipping MALDIquant benchmark for $(case.name) to prevent SIGKILL OOM.")
                push!(results, (
                    Dataset = case.name,
                    FileSize_MB = file_size,
                    Library = "MALDIquant",
                    LoadTime_s = NaN,
                    LoadMem_MB = NaN,
                    AvgSliceTime_s = NaN,
                    SliceStd_s = 0.0,
                    Status = "OOM / Skipped"
                ))
            else
                println("[R - MALDIquant (centroid=TRUE)]")
                temp_r_csv = joinpath(RESULTS_DIR, "temp_r_mq.csv")
                r_script_path = joinpath(@__DIR__, "benchmark_r.R")
                
                try
                    cmd = `Rscript $r_script_path maldiquant $(case.filepath) $(case.mz_value) $(case.mz_tolerance) $NUM_REPETITIONS $temp_r_csv`
                    process_status = run(cmd, wait=true)
                    
                    if success(process_status) && isfile(temp_r_csv)
                        r_df = CSV.read(temp_r_csv, DataFrame)
                        push!(results, (
                            Dataset = case.name,
                            FileSize_MB = file_size,
                            Library = "MALDIquant",
                            LoadTime_s = r_df.LoadTime_s[1],
                            LoadMem_MB = r_df.LoadMem_MB[1],
                            AvgSliceTime_s = r_df.AvgSliceTime_s[1],
                            SliceStd_s = 0.0, # R script simplified to mean
                            Status = "Success"
                        ))
                        rm(temp_r_csv, force=true)
                    else
                        throw(ErrorException("R process failed"))
                    end
                catch e
                    @warn "MALDIquant benchmark failed/OOMed: $e"
                    push!(results, (
                        Dataset = case.name,
                        FileSize_MB = file_size,
                        Library = "MALDIquant",
                        LoadTime_s = NaN,
                        LoadMem_MB = NaN,
                        AvgSliceTime_s = NaN,
                        SliceStd_s = 0.0,
                        Status = "OOM / Failed"
                    ))
                end
            end

            # ------------------------------------------------------------
            # 3. Run R Cardinal (Isolated process)
            # ------------------------------------------------------------
            println("[R - Cardinal (matter lazy-loading)]")
            temp_r_csv_card = joinpath(RESULTS_DIR, "temp_r_card.csv")
            r_script_path = joinpath(@__DIR__, "benchmark_r.R")
            
            try
                cmd = `Rscript $r_script_path cardinal $(case.filepath) $(case.mz_value) $(case.mz_tolerance) $NUM_REPETITIONS $temp_r_csv_card`
                process_status = run(cmd, wait=true)
                
                if success(process_status) && isfile(temp_r_csv_card)
                    r_df = CSV.read(temp_r_csv_card, DataFrame)
                    push!(results, (
                        Dataset = case.name,
                        FileSize_MB = file_size,
                        Library = "Cardinal",
                        LoadTime_s = r_df.LoadTime_s[1],
                        LoadMem_MB = r_df.LoadMem_MB[1],
                        AvgSliceTime_s = r_df.AvgSliceTime_s[1],
                        SliceStd_s = 0.0,
                        Status = "Success"
                    ))
                    rm(temp_r_csv_card, force=true)
                else
                    throw(ErrorException("R process failed"))
                end
            catch e
                @warn "Cardinal benchmark failed/OOMed: $e"
                push!(results, (
                    Dataset = case.name,
                    FileSize_MB = file_size,
                    Library = "Cardinal",
                    LoadTime_s = NaN,
                    LoadMem_MB = NaN,
                    AvgSliceTime_s = NaN,
                    SliceStd_s = 0.0,
                    Status = "OOM / Failed"
                ))
            end
        end
    end

    # Calculate Derived Academic Metrics
    # Throughput, MemEfficiency, FoldSpeedup (rel. to Cardinal), LoadBandwidth
    results.Throughput_SPS = 1.0 ./ results.AvgSliceTime_s
    results.MemEff_MB_GB = results.LoadMem_MB ./ (results.FileSize_MB ./ 1024.0)
    results.LoadBandwidth_MBs = results.FileSize_MB ./ results.LoadTime_s

    # Fold Speedup calculation using Cardinal as Baseline
    results.Fold_Speedup = fill(NaN, nrow(results))
    for ds in unique(results.Dataset)
        ds_rows = results[results.Dataset .== ds, :]
        cardinal_row = ds_rows[ds_rows.Library .== "Cardinal", :]
        julia_row = ds_rows[ds_rows.Library .== "JuliaMSI", :]
        if nrow(cardinal_row) > 0 && nrow(julia_row) > 0 && cardinal_row.Status[1] == "Success"
            card_time = cardinal_row.AvgSliceTime_s[1]
            julia_time = julia_row.AvgSliceTime_s[1]
            # Speedup = Cardinal Time / JuliaMSI Time
            speedup = card_time / julia_time
            results[ (results.Dataset .== ds) .& (results.Library .== "JuliaMSI"), :Fold_Speedup ] .= speedup
        end
    end

    csv_path = joinpath(RESULTS_DIR, "cross_language_benchmark.csv")
    CSV.write(csv_path, results)
    println("\nDetailed cross-language results saved to $csv_path")
    
    if nrow(results) > 0
        plot_cross_results(results)
    end
    
    return results
end

function plot_cross_results(df::DataFrame)
    unique_datasets = unique(df.Dataset)
    if isempty(unique_datasets) return end

    fig = Figure(size=(1800, 1200), fontsize=22)
    
    libraries = ["JuliaMSI", "MALDIquant", "Cardinal"]
    colors = ["#10b981", "#3b82f6", "#8b5cf6"] # green, blue, purple
    x_pos = 1:length(unique_datasets)
    
    # ---------- Panel A: Load Time (Log Scale) ----------
    ax1 = Axis(fig[1, 1],
               title="Panel A: Metadata/File Load Time",
               ylabel="Time (seconds) - Log Scale",
               yscale=log10,
               xticks=(x_pos, unique_datasets), xticklabelrotation=π/8)
    
    # ---------- Panel B: Memory Consumption (Log Scale) ----------
    ax2 = Axis(fig[1, 2],
               title="Panel B: Max RAM Footprint",
               ylabel="Memory (MB) - Log Scale",
               yscale=log10,
               xticks=(x_pos, unique_datasets), xticklabelrotation=π/8)
               
    # ---------- Panel C: Slice Generation Time (Log Scale) ----------
    ax3 = Axis(fig[2, 1:2],
               title="Panel C: Slice Generation Time (Mean ± SD)",
               ylabel="Time (seconds) - Log Scale",
               yscale=log10,
               xticks=(x_pos, unique_datasets), xticklabelrotation=π/8)

    width = 0.25
    offsets = [-width, 0.0, width]
    
    for (idx, lib) in enumerate(libraries)
        lib_df = df[df.Library .== lib, :]
        if isempty(lib_df) continue end
        
        load_times = Float64[]
        load_mems = Float64[]
        slice_means = Float64[]
        slice_stds = Float64[]
        
        for ds in unique_datasets
            row = lib_df[lib_df.Dataset .== ds, :]
            if nrow(row) > 0 && row.Status[1] == "Success"
                push!(load_times, coalesce(row.LoadTime_s[1], 1e-3))
                push!(load_mems, coalesce(row.LoadMem_MB[1], 1e-1))
                push!(slice_means, coalesce(row.AvgSliceTime_s[1], 1e-4))
                push!(slice_stds, coalesce(row.SliceStd_s[1], 0.0))
            else
                # Default tiny values for log plot positioning of missing runs
                push!(load_times, 1e-3)
                push!(load_mems, 1e-1)
                push!(slice_means, 1e-4)
                push!(slice_stds, 0.0)
            end
        end
        
        barplot!(ax1, x_pos .+ offsets[idx], load_times, color=colors[idx], width=width, label=lib)
        barplot!(ax2, x_pos .+ offsets[idx], load_mems, color=colors[idx], width=width, label=lib)
        
        barplot!(ax3, x_pos .+ offsets[idx], slice_means, color=colors[idx], width=width, label=lib)
        errorbars!(ax3, x_pos .+ offsets[idx], slice_means, slice_stds, color=:black, whiskerwidth=10)
    end
    
    axislegend(ax1, position=:lt)
    axislegend(ax2, position=:lt)
    axislegend(ax3, position=:rt)
    
    save(joinpath(RESULTS_DIR, "cross_language_comparison.png"), fig)
    println("Comparison dashboard saved to $(joinpath(RESULTS_DIR, "cross_language_comparison.png"))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cross_benchmarks()
end
