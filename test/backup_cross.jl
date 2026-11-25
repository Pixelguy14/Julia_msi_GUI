# test/benchmark_comparison.jl

# ===================================================================
# Cross-Language Performance Benchmark Suite: JuliaMSI vs R (MALDIquant & Cardinal)
# ===================================================================
# This script runs the performance benchmarks for:
# 1. JuliaMSI (New Engine)
# 2. R - MALDIquant (isolated run)
# 3. R - Cardinal (isolated run)
#
# It saves the results to `test/results/cross_language_benchmark.csv` and
# creates a comparison plot.
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
const NUM_REPETITIONS = 10 # Number of slice repetitions for averaging

const BENCHMARK_CASES = [
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_LA-ESI/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML", 436.04, 0.1, "Thaliana Leaf"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Chilli/ltpmsi-chilli.imzML", 420.0, 0.1, "Chilli Pepper"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_DESI/40TopL,10TopR,30BottomL,20BottomR/40TopL,10TopR,30BottomL,20BottomR-centroid.imzML", 885.5, 0.1, "Human Colon Cancer"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML", 716.053, 0.1, "Mouse Urinary Bladder"),
    #BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Leafs/CE1_Leaf_R3.imzML", 306.1, 0.1, "Leaf"), # Este leaf es diferente de taliana!
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

    # Check if Rscript is available
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
        mean_time_new_s = mean(b_slice_new.times) / 1e9

        close(msi_data)
        flush_memory()

        # Record JuliaMSI results
        push!(results, (
            Dataset = case.name,
            FileSize_MB = file_size,
            Library = "JuliaMSI",
            LoadTime_s = load_time_new_s,
            LoadMem_MB = mem_load_new_mb,
            AvgSliceTime_s = mean_time_new_s,
            Status = "Success"
            ))

        # ------------------------------------------------------------
        # 2. Run R MALDIquant (Isolated process)
        # ------------------------------------------------------------
        if has_r
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
                        Status = "Success"
                        ))
                    rm(temp_r_csv, force=true)
                else
                    throw(ErrorException("R process failed with exit code $(process_status.exitcode)"))
                end
                catch e
                @warn "MALDIquant benchmark failed/OOMed for $(case.name): $e"
                    push!(results, (
                        Dataset = case.name,
                        FileSize_MB = file_size,
                        Library = "MALDIquant",
                        LoadTime_s = NaN,
                        LoadMem_MB = NaN,
                        AvgSliceTime_s = NaN,
                        Status = "OOM / Failed"
                        ))
                end

                # ------------------------------------------------------------
                # 3. Run R Cardinal (Isolated process)
                # ------------------------------------------------------------
                println("[R - Cardinal (matter lazy-loading)]")
                temp_r_csv_card = joinpath(RESULTS_DIR, "temp_r_card.csv")

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
                            Status = "Success"
                            ))
                        rm(temp_r_csv_card, force=true)
                    else
                        throw(ErrorException("R process failed with exit code $(process_status.exitcode)"))
                    end
                    catch e
                    @warn "Cardinal benchmark failed/OOMed for $(case.name): $e"
                        push!(results, (
                            Dataset = case.name,
                            FileSize_MB = file_size,
                            Library = "Cardinal",
                            LoadTime_s = NaN,
                            LoadMem_MB = NaN,
                            AvgSliceTime_s = NaN,
                            Status = "OOM / Failed"
                            ))
                    end
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
            # Filter datasets
            unique_datasets = unique(df.Dataset)
            if isempty(unique_datasets) return end

            fig = Figure(size=(1800, 1200), fontsize=20)

            libraries = ["JuliaMSI", "MALDIquant", "Cardinal"]
            colors = ["#10b981", "#3b82f6", "#8b5cf6"] # green, blue, purple

            x_pos = 1:length(unique_datasets)

            # Axis 1: Slice Generation Time
            ax1 = Axis(fig[1, 1],
                       title="Average Slice Generation Time (Lower is Better)",
                       ylabel="Time (seconds)",
                       xticks=(x_pos, unique_datasets), xticklabelrotation=π/8)

            # Axis 2: RAM overhead during load
            ax2 = Axis(fig[1, 2],
                       title="Memory Consumption during Load (Lower is Better)",
                       ylabel="Memory (MB)",
                       xticks=(x_pos, unique_datasets), xticklabelrotation=π/8)

            # Plot grouped bars
            width = 0.25
            offsets = [-width, 0.0, width]

            for (idx, lib) in enumerate(libraries)
                lib_df = df[df.Library .== lib, :]
                if isempty(lib_df) continue end

                # Match rows with unique_datasets order
                slice_times = Float64[]
                load_mems = Float64[]

                for ds in unique_datasets
                    row = lib_df[lib_df.Dataset .== ds, :]
                    if nrow(row) > 0 && row.Status[1] == "Success"
                        push!(slice_times, coalesce(row.AvgSliceTime_s[1], 0.0))
                        push!(load_mems, coalesce(row.LoadMem_MB[1], 0.0))
                    else
                        push!(slice_times, 0.0)
                        push!(load_mems, 0.0)
                    end
                end

                barplot!(ax1, x_pos .+ offsets[idx], slice_times, color=colors[idx], width=width, label=lib)
                barplot!(ax2, x_pos .+ offsets[idx], load_mems, color=colors[idx], width=width, label=lib)
            end

            axislegend(ax1, position=:rt)
            axislegend(ax2, position=:rt)

            save(joinpath(RESULTS_DIR, "cross_language_comparison.png"), fig)
            println("Comparison dashboard saved to $(joinpath(RESULTS_DIR, "cross_language_comparison.png"))")
        end

        if abspath(PROGRAM_FILE) == @__FILE__
            run_cross_benchmarks()
        end
