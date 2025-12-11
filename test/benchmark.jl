# test/benchmark.jl

# ===================================================================
# Performance Benchmark Suite for JuliaMSI
# ===================================================================
# This script compares the performance of the new `JuliaMSI` pipeline
# against the old `julia_mzML_imzML` library.
#
# It measures the time and memory allocation for loading an `.imzML`
# file and generating an image slice multiple times.
#
# Instructions:
# 1. Fill in the placeholder paths in the "CONFIG" section below.
# 2. Run the script from the project's root directory:
#    julia --project=. test/benchmark.jl
# 3. Check `test/results/` for `benchmark_results.csv` and plots.
# ===================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) # Activate main project environment

using DataFrames
using CSV
using CairoMakie
using ColorSchemes
using Statistics
using Printf
using JSON
using MSI_src # Import the new library
using .MSI_src: get_mz_slice, quantize_intensity, save_bitmap, get_multiple_mz_slices
using julia_mzML_imzML

# Remember, to run this, you must enter pkg mode in julia and add this library:
# add https://github.com/CINVESTAV-LABI/julia_mzML_imzML

# ===================================================================
# CONFIG: PLEASE FILL IN YOUR FILE PATHS AND PARAMETERS
# ===================================================================

struct BenchmarkCase
    filepath::String
    mz_value::Float64
    mz_tolerance::Float64
    name::String # New field for a descriptive name
end

const BENCHMARK_CASES = [
    # Add your benchmark cases here. Example:
    # BenchmarkCase("/path/to/your/small_file.imzML", 309.06, 0.1),
    # BenchmarkCase("/path/to/your/medium_file.imzML", 896.0, 1.0),
    # BenchmarkCase("/path/to/your/large_file.imzML", 100.0, 0.1),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_LA-ESI/180817_NEG_Thaliana_Leaf_bottom_1_0841.imzML",116.07,0.1, "Thaliana Leaf"),
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_LTP/ltpmsi-chilli.imzML",420,0.1, "Chilli Pepper"), # Chilli
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_DESI/ColAd_Individual/40TopL,10TopR,30BottomL,20BottomR/40TopL,10TopR,30BottomL,20BottomR-centroid.imzML",885.5,0.1, "Colon Cancer Human"), # Human Cancer
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/imzML_AP_SMALDI/HR2MSImouseurinarybladderS096.imzML", 716.053,0.1, "Mouse Urinary Bladder"), # Mouse bladder
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/Liv2_imzML_TIMSConvert-selected/Liv2.imzML",796.18,0.1, "Liver Cut"), #Lib2
    BenchmarkCase("/home/pixel/Documents/Cinvestav_2025/Analisis/salida/Stomach_DHB_uncompressed.imzML",804.3,0.1, "Mouse Stomach"), # Mouse Stomach
]

const NUM_REPETITIONS = 50 # Number of times to generate the image for averaging
const BULK_MZ_VARIATIONS = 50 # Number of m/z variations to generate for the bulk benchmark (e.g., mz, mz+0.1, ..., mz+(N-1)*0.1)
const RESULTS_DIR = joinpath(@__DIR__, "results")

# Helper function to generate m/z variations for bulk processing
function get_mz_variation_vector(base_mz::Float64, num_variations::Int)
    return [base_mz + i * 0.1 for i in 0:(num_variations - 1)]
end

# Helper function to get the combined size of .imzML and .ibd files
function get_total_file_size_mb(filepath::String)
    imzml_size = isfile(filepath) ? filesize(filepath) : 0
    ibd_path = replace(filepath, r"\.(imzML|imzml)$" => ".ibd")
    ibd_size = isfile(ibd_path) ? filesize(ibd_path) : 0
    return round((imzml_size + ibd_size) / 1024^2, digits=2)
end

# ===================================================================
# BENCHMARK: JuliaMSI (New Library)
# ===================================================================

"""
    benchmark_julianew(case::BenchmarkCase)

Runs the benchmark for the modern JuliaMSI library.
"""
function benchmark_julianew(case::BenchmarkCase)
    println("--- Benchmarking JuliaMSI (New) on $(basename(case.filepath)) ---")
    
    # 1. Load the data
    load_stats = @timed OpenMSIData(case.filepath)
    msi_data = load_stats.value
    load_time = load_stats.time
    load_mem = load_stats.bytes / 1024^2 # Convert to MB

    image_times = []
    image_mems = []
    local final_image # Declare final_image outside the loop

    # 2. Generate image multiple times
    for i in 1:NUM_REPETITIONS
        print("\rRepetition $i/$NUM_REPETITIONS...")
        image_stats = @timed get_mz_slice(msi_data, case.mz_value, case.mz_tolerance)
        final_image = image_stats.value
        push!(image_times, image_stats.time)
        push!(image_mems, image_stats.bytes / 1024^2)
    end
    println("\nDone.")

    # 3. Save one resulting image for validation
    output_path_bitmap = joinpath(RESULTS_DIR, "benchmark_JuliaMSI_$(basename(case.filepath)).bmp")
    # Quantize and save the image
    quantized_image, _ = quantize_intensity(final_image, 20) # 256 levels
    save_bitmap(output_path_bitmap, quantized_image, MSI_src.ViridisPalette)
    println("Saved validation image to $output_path_bitmap")

    # 4. Calculate statistics
    println("DEBUG (New): NUM_REPETITIONS = $NUM_REPETITIONS")
    println("DEBUG (New): load_time = $load_time")
    println("DEBUG (New): image_times = $(image_times)")
    println("DEBUG (New): length(image_times) = $(length(image_times))")
    avg_image_time = mean(image_times)
    avg_image_mem = mean(image_mems) # Added missing calculation
    println("DEBUG (New): avg_image_time = $avg_image_time")
    println("DEBUG (New): sum(image_times) = $(sum(image_times))")
    total_time = load_time + sum(image_times)
    println("DEBUG (New): calculated total_time = $total_time")

    return (
        library="JuliaMSI",
        filepath=basename(case.filepath),
        file_size_mb=get_total_file_size_mb(case.filepath), # Use helper function
        load_time_s=load_time,
        avg_image_time_s=avg_image_time,
        total_time_s=total_time,
        load_memory_mb=load_mem,
        avg_image_memory_mb=avg_image_mem,
        name=case.name # Added case name to results
    )
end

"""
    benchmark_julianew_bulk(case::BenchmarkCase)

Runs the benchmark for the modern JuliaMSI library's bulk slice generation.
"""
function benchmark_julianew_bulk(case::BenchmarkCase)
    println("--- Benchmarking JuliaMSI (Bulk) on $(basename(case.filepath)) ---")

    # Generate multiple m/z values for bulk processing
    mz_values_for_bulk = get_mz_variation_vector(case.mz_value, BULK_MZ_VARIATIONS)
    println("  Generating $(length(mz_values_for_bulk)) m/z variations: $(mz_values_for_bulk)")
    
    # 1. Load the data
    load_stats = @timed OpenMSIData(case.filepath)
    msi_data = load_stats.value
    load_time = load_stats.time
    load_mem = load_stats.bytes / 1024^2 # Convert to MB

    image_times = []
    image_mems = []
    local final_slices_dict # Declare outside the loop

    # 2. Generate multiple images in bulk (repeated for averaging)
    for i in 1:1
        print("\rRepetition $i/$NUM_REPETITIONS for bulk...")
        # Assume MSI_src.get_multiple_mz_slices exists and returns Dict{Real, Matrix{Float64}}
        image_stats = @timed MSI_src.get_multiple_mz_slices(msi_data, mz_values_for_bulk, case.mz_tolerance)
        final_slices_dict = image_stats.value
        push!(image_times, image_stats.time)
        
        # Calculate memory for all slices combined
        total_slices_mem = sum(Base.summarysize(slice) for slice in values(final_slices_dict)) / 1024^2
        push!(image_mems, total_slices_mem)
    end
    println("\nDone.")

    # No validation image saving for bulk due to multiple slices

    # 3. Calculate statistics
    avg_image_time = mean(image_times)
    avg_image_mem = mean(image_mems)
    total_time = load_time + sum(image_times) # Total time includes all repetitions

    return (
        library="JuliaMSIBulk",
        filepath=basename(case.filepath),
        file_size_mb=get_total_file_size_mb(case.filepath), # Use helper function
        load_time_s=load_time,
        avg_image_time_s=avg_image_time,
        total_time_s=total_time,
        load_memory_mb=load_mem,
        avg_image_memory_mb=avg_image_mem,
        name=case.name # Added case name to results
    )
end


"""
    benchmark_juliaold(case::BenchmarkCase)

Runs the benchmark for the old library using dynamic environment switching.
"""
function benchmark_juliaold(case::BenchmarkCase)
    println("--- Benchmarking julia_mzML_imzML (Old) on $(basename(case.filepath)) ---") 
    try
        # --- 1. Load Data ---
        load_stats = @timed LoadImzml(case.filepath)
        msi_data = load_stats.value
        load_time = load_stats.time
        load_mem = load_stats.bytes / 1024^2 # Convert to MB

        image_times = []
        image_mems = []
        local last_image_stats_value # Declare outside the loop

        # --- 2. Generate Image Repeatedly ---
        for i in 1:NUM_REPETITIONS
            print("\rRepetition $i/$NUM_REPETITIONS...")
            image_stats = @timed GetSlice(msi_data, case.mz_value, case.mz_tolerance)
            last_image_stats_value = image_stats.value # Store the value
            push!(image_times, image_stats.time)
            push!(image_mems, image_stats.bytes / 1024^2)
        end
        println("\nDone.")

        # --- 3. Save Validation Image ---
        final_image = last_image_stats_value
        output_path_bitmap = joinpath(RESULTS_DIR, "benchmark_old_lib_$(basename(case.filepath)).bmp")
        SaveBitmap(output_path_bitmap, IntQuant(final_image), ViridisPalette)
        println("Saved validation image to $output_path_bitmap")

            # --- 4. Calculate Averages and Totals ---
            println("DEBUG (Old): NUM_REPETITIONS = $NUM_REPETITIONS")
            println("DEBUG (Old): load_time = $load_time")
            println("DEBUG (Old): image_times = $(image_times)")
            println("DEBUG (Old): length(image_times) = $(length(image_times))")
            avg_image_time = mean(image_times)
            println("DEBUG (Old): avg_image_time = $avg_image_time")
            println("DEBUG (Old): sum(image_times) = $(sum(image_times))")
            total_time = load_time + sum(image_times)
            println("DEBUG (Old): calculated total_time = $total_time")
            avg_image_mem = mean(image_mems) # Added missing calculation
        return (
            library="julia_mzML_imzML",
            filepath=basename(case.filepath),
            file_size_mb=get_total_file_size_mb(case.filepath), # Use helper function
            load_time_s=load_time,
            avg_image_time_s=avg_image_time,
            total_time_s=total_time,
        load_memory_mb=load_mem,
        avg_image_memory_mb=avg_image_mem,
        name=case.name # Added case name to results
    )
        
    catch e
        println("ERROR benchmarking old library: $e")
        println("Make sure julia_mzML_imzML is installed and OLD_LIB_ENV_PATH is set correctly.")
        return nothing
    end
end

function benchmark_multi_file_walkthrough(cases::Vector{BenchmarkCase})
    println("--- Multi-File Walkthrough Benchmark ---")
    println("Simulating workflow: load file -> generate 1 image -> close -> next file")
    
    results = []
    
    for case in cases
        if !isfile(case.filepath)
            @warn "File not found: $(case.filepath). Skipping."
            continue
        end
        
        println("\nProcessing: $(case.name)")
        
        # --- NEW LIBRARY (JuliaMSI) ---
        println("  [JuliaMSI]")
        # Time includes metadata loading + single image generation
        new_total_stats = @timed begin
            data = OpenMSIData(case.filepath)
            img = get_mz_slice(data, case.mz_value, case.mz_tolerance)
            close(data) # Explicitly release resources
        end
        
        # --- OLD LIBRARY (julia_mzML_imzML) ---
        println("  [julia_mzML_imzML]")
        old_total_stats = try
            @timed begin
                data = LoadImzml(case.filepath) # Loads ENTIRE dataset into RAM
                img = GetSlice(data, case.mz_value, case.mz_tolerance)
                # Note: No explicit close function in old library
            end
        catch e
            @warn "Old library failed on $(basename(case.filepath)): $e"
            (time=Inf, bytes=Inf, value=nothing)
        end
        
        push!(results, (
            name=case.name,
            file_size_mb=get_total_file_size_mb(case.filepath),
            julianew_time_s=new_total_stats.time,
            julianew_memory_mb=new_total_stats.bytes/1024^2,
            juliaold_time_s=old_total_stats.time,
            juliaold_memory_mb=old_total_stats.bytes/1024^2,
            speedup=old_total_stats.time / new_total_stats.time
        ))
        
        # Force garbage collection between files to level the playing field
        GC.gc()
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0)
        end
    end
    
    return results
end

function benchmark_real_world_pipeline(cases::Vector{BenchmarkCase})
    println("--- REAL-WORLD PIPELINE BENCHMARK ---")
    println("Simulating: Load → Multiple Slices → Process → Save → Next File")
    
    # Realistic m/z values (similar to your UI)
    mz_values = [116.07, 420.0, 885.5, 716.053, 796.18, 804.3]
    tolerance = 0.1
    results = []
    
    for case in cases
        if !isfile(case.filepath)
            @warn "File not found: $(case.filepath). Skipping."
            continue
        end
        
        println("\nProcessing: $(case.name)")
        
        # --- NEW LIBRARY (JuliaMSI) with bulk processing ---
        println("  [JuliaMSI - Bulk Processing]")
        new_stats = try
            @timed begin
                # 1. Load (lazy, metadata only)
                data = OpenMSIData(case.filepath)
                
                # 2. Generate ALL slices in one bulk operation
                slice_dict = get_multiple_mz_slices(data, mz_values, tolerance)
                
                # 3. Process each slice (simulate quantization, saving, etc.)
                for (mass, slice) in slice_dict
                    # Simulate minimal processing for fair comparison
                    # In real pipeline: TrIQ, median filter, save bitmap
                    processed = quantize_intensity(slice, 256)[1]
                end
                
                # 4. Close and cleanup
                close(data)
            end
        catch e
            @warn "New library failed: $e"
            (time=Inf, bytes=Inf, value=nothing)
        end
        
        # --- OLD LIBRARY (sequential processing) ---
        println("  [Old Library - Sequential]")
        old_stats = try
            @timed begin
                # 1. Load (everything into RAM)
                data = LoadImzml(case.filepath)
                
                # 2. Generate slices one by one
                for mass in mz_values
                    slice = GetSlice(data, mass, tolerance)
                    # Same processing simulation
                    processed = IntQuant(slice)
                end
                
                # 3. Old library doesn't have explicit close
                # Data stays in RAM until GC
            end
        catch e
            @warn "Old library failed: $e"
            (time=Inf, bytes=Inf, value=nothing)
        end
        
        # Memory cleanup between files (mirroring your pipeline)
        GC.gc(true)
        if Sys.islinux()
            ccall(:malloc_trim, Int32, (Int32,), 0)
        end
        
        push!(results, (
            name=case.name,
            file_size_mb=get_total_file_size_mb(case.filepath),
            n_slices=length(mz_values),
            julianew_time_s=new_stats.time,
            julianew_memory_mb=new_stats.bytes/1024^2,
            juliaold_time_s=old_stats.time,
            juliaold_memory_mb=old_stats.bytes/1024^2,
            speedup=old_stats.time / new_stats.time,
            # NEW METRICS for your pipeline:
            memory_reduction_pct=100 * (old_stats.bytes - new_stats.bytes) / old_stats.bytes,
            slices_per_second_new=length(mz_values) / new_stats.time,
            slices_per_second_old=length(mz_values) / old_stats.time
        ))
    end
    
    return results
end

# ===================================================================
# Main Runner
function run_all_benchmarks(skip_benchmark::Bool=false)
    results_csv_path = joinpath(RESULTS_DIR, "benchmark_results.csv")
    multi_file_walkthrough_csv_path = joinpath(RESULTS_DIR, "multi_file_walkthrough.csv")
    real_world_pipeline_csv_path = joinpath(RESULTS_DIR, "real_world_pipeline.csv")

    local all_results::DataFrame
    local walkthrough_df::DataFrame = DataFrame() # Initialize as empty
    local pipeline_df::DataFrame = DataFrame() # Initialize as empty

    if skip_benchmark && isfile(results_csv_path) && isfile(multi_file_walkthrough_csv_path) && isfile(real_world_pipeline_csv_path)
        @info "Skipping benchmark run. Loading results from $results_csv_path"
        all_results = CSV.read(results_csv_path, DataFrame)
        walkthrough_df = CSV.read(multi_file_walkthrough_csv_path, DataFrame)
        pipeline_df = CSV.read(real_world_pipeline_csv_path, DataFrame)
    else
        if isempty(BENCHMARK_CASES)
            @warn "No benchmark cases found. Please fill in the `BENCHMARK_CASES` constant in `test/benchmark.jl`."
            return
        end
        
        all_results = DataFrame()

        for case in BENCHMARK_CASES
            if !isfile(case.filepath)
                @warn "File not found: $(case.filepath). Skipping."
                continue
            end

            # Run for new library
            new_results = benchmark_julianew(case)
            push!(all_results, new_results, cols=:union)

            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end

            # Run for new library (Bulk)
            new_bulk_results = benchmark_julianew_bulk(case)
            push!(all_results, new_bulk_results, cols=:union)

            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end

            # Run for old library
            old_results = benchmark_juliaold(case)
            if old_results !== nothing
                push!(all_results, old_results, cols=:union)
            end
            
            GC.gc() # Trigger garbage collection
            if Sys.islinux()
                ccall(:malloc_trim, Int32, (Int32,), 0) # Ensure julia returns the freed memory to OS
            end
        end

        # Save results to CSV (only if benchmark was actually run)
        CSV.write(results_csv_path, all_results)
        println("\nBenchmark results saved to $results_csv_path")

        # --- Multi-File Walkthrough Benchmark ---
        println("\n" * "="^60)
        println("MULTI-FILE WALKTHROUGH BENCHMARK")
        println("="^60)

        walkthrough_results = benchmark_multi_file_walkthrough(BENCHMARK_CASES)

        # Create a summary DataFrame
        walkthrough_df = DataFrame(
            name=[r.name for r in walkthrough_results],
            file_size_mb=[r.file_size_mb for r in walkthrough_results],
            JuliaMSI_time_s=[r.julianew_time_s for r in walkthrough_results],
            Old_time_s=[r.juliaold_time_s for r in walkthrough_results],
            speedup=[r.speedup for r in walkthrough_results],
            JuliaMSI_mem_mb=[r.julianew_memory_mb for r in walkthrough_results],
            Old_mem_mb=[r.juliaold_memory_mb for r in walkthrough_results]
        )

        println("\nMulti-file walkthrough results:")
        show(walkthrough_df, allrows=true)

        # Save to CSV
        CSV.write(multi_file_walkthrough_csv_path, walkthrough_df)

        # --- REAL-WORLD PIPELINE BENCHMARK ---
        println("\n" * "="^60)
        println("REAL-WORLD PIPELINE BENCHMARK")
        println("="^60)

        pipeline_results = benchmark_real_world_pipeline(BENCHMARK_CASES)

        # Create a summary DataFrame
        pipeline_df = DataFrame(
            name=[r.name for r in pipeline_results],
            file_size_mb=[r.file_size_mb for r in pipeline_results],
            n_slices=[r.n_slices for r in pipeline_results],
            JuliaMSI_time_s=[r.julianew_time_s for r in pipeline_results],
            Old_time_s=[r.juliaold_time_s for r in pipeline_results],
            speedup=[r.speedup for r in pipeline_results],
            JuliaMSI_mem_mb=[r.julianew_memory_mb for r in pipeline_results],
            Old_mem_mb=[r.juliaold_memory_mb for r in pipeline_results],
            memory_reduction_pct=[r.memory_reduction_pct for r in pipeline_results],
            slices_per_second_new=[r.slices_per_second_new for r in pipeline_results],
            slices_per_second_old=[r.slices_per_second_old for r in pipeline_results]
        )

        println("\nReal-world pipeline results:")
        show(pipeline_df, allrows=true)

        # Save to CSV
        CSV.write(real_world_pipeline_csv_path, pipeline_df)
    end

    # Generate and save plots (always happens)
    generate_plots(all_results)
    plot_walkthrough_results(walkthrough_df)
    plot_real_world_pipeline_results(pipeline_df)
end

# ===================================================================
# PLOTTING - Professional comparative visualizations
# ===================================================================    
function generate_plots(df::DataFrame)
    if isempty(df)
        println("Cannot generate plots from empty dataframe.")
        return
    end

    unique_names = unique(df.name) # Use the new name field
    num_files = length(unique_names)
    
    # Calculate positions for dodged bars
    positions = 1:num_files # Define positions here
    
    # Generate A-Z labels for X-axis
    alphabet_labels = [string(Char('A' + i - 1)) for i in 1:num_files]
    
    # Use a scientific color palette
    colors = Dict(
        "JuliaMSI" => "#2E86AB",  # Professional blue
        "julia_mzML_imzML" => "#A23B72",  # Professional magenta
        "JuliaMSIBulk" => "#2CA02C" # Professional green for bulk
    )
    
    # Shorter labels for display
    library_labels = Dict(
        "JuliaMSI" => "New",
        "julia_mzML_imzML" => "Old",
        "JuliaMSIBulk" => "New Bulk"
    )
    
    # Define specific colors for the speedup metrics for each library
    speedup_colors_new = ["#1E88E5", "#43A047", "#FFB300", "#E53935"] # Blue, Green, Yellow, Red for JuliaMSI
    speedup_colors_bulk = ["#0056B3", "#008000", "#FFA500", "#CC0000"] # Slightly darker/different hues for JuliaMSIBulk
    
    # Set modern theme
    set_theme!(Theme(
        fontsize = 72, # General font size (1.5x of 48)
        font = "Helvetica", # Changed to Helvetica
        Axis = (
            backgroundcolor = :white,
            spinewidth = 1.5,
            xgridvisible = false,
            ygridvisible = true,
            ygridcolor = (:gray, 0.2),
            ygridwidth = 0.5,
            xticklabelrotation = π/8,
            xticklabalign = (:right, :center),
            yticklabalign = (:right, :center),
            xlabelsize = 72, # Axis label size (1.5x of 48)
            ylabelsize = 72, # Axis label size (1.5x of 48)
            titlesize = 72, # Axis title size (1.5x of 48)
            xticklabelsize = 60, # Tick label size (1.5x of 40)
            yticklabelsize = 60, # Tick label size (1.5x of 40)
            titlefont = :bold, # Ensured bold
        ),
        Legend = (
            backgroundcolor = (:white, 0.95),
            framecolor = :gray,
            framewidth = 1,
            padding = (10, 10, 5, 5),
            labelsize = 72, # Legend label size (1.5x of 48)
            titlesize = 72, # Legend title size (1.5x of 48)
        )
    ))

    # ===================================================================
    # Plot 1: Main Performance Dashboard
    # ===================================================================
    fig1 = Figure(size = (3600, 2400)) # Explicitly set figure size
    
    # Define grid layout with specific areas
    g1 = fig1[2, 1:2] = GridLayout()
    
    # Total Time
    ax1 = Axis(g1[1, 1],
        title = "Total Time",
        ylabel = "Seconds", # Changed to linear scale, removed "(log scale)"
        xticks = (positions, alphabet_labels) # Use alphabetical labels for x-axis
    )
    
    # Average Image Time
    ax2 = Axis(g1[1, 2],
        title = "Image Generation Time",
        ylabel = "Seconds",
        yticks = LinearTicks(5),
        xticks = (positions, alphabet_labels) # Use alphabetical labels for x-axis
    )
    
    # Loading Memory
    ax3 = Axis(g1[2, 1],
        title = "Loading Memory",
        ylabel = "MB",
        yticks = LinearTicks(5),
        xticks = (positions, alphabet_labels) # Use alphabetical labels for x-axis
    )
    
    # Image Generation Memory
    ax4 = Axis(g1[2, 2],
        title = "Image Generation Memory",
        ylabel = "MB",
        yticks = LinearTicks(5),
        xticks = (positions, alphabet_labels) # Use alphabetical labels for x-axis
    )

    bar_width = 0.25 # Reduced bar width
    x_offsets = [-0.3, 0, 0.3] # Offsets for 3 bars to avoid overlap
    
    # Arrays to store speedup factors
    total_time_speedup_new = zeros(num_files)
    image_time_speedup_new = zeros(num_files)
    load_mem_speedup_new = zeros(num_files)
    img_mem_speedup_new = zeros(num_files)

    total_time_speedup_bulk = zeros(num_files)
    image_time_speedup_bulk = zeros(num_files)
    load_mem_speedup_bulk = zeros(num_files)
    img_mem_speedup_bulk = zeros(num_files)
    
    # Initialize a list to hold rows for the summary DataFrame
    summary_df_rows = []

    for (i, name_for_df_filter) in enumerate(unique_names) # Iterate through unique names
        # Get data for this file
        new_data = df[(df.name .== name_for_df_filter) .& (df.library .== "JuliaMSI"), :]
        bulk_data = df[(df.name .== name_for_df_filter) .& (df.library .== "JuliaMSIBulk"), :]
        old_data = df[(df.name .== name_for_df_filter) .& (df.library .== "julia_mzML_imzML"), :]
        
        # Plot total time
        barplot!(ax1, [i + x_offsets[1], i + x_offsets[2], i + x_offsets[3]],
            [new_data.total_time_s[1], bulk_data.total_time_s[1], old_data.total_time_s[1]],
            color = [colors["JuliaMSI"], colors["JuliaMSIBulk"], colors["julia_mzML_imzML"]],
            width = bar_width
        )
        
        # Plot average image time
        barplot!(ax2, [i + x_offsets[1], i + x_offsets[2], i + x_offsets[3]],
            [new_data.avg_image_time_s[1], bulk_data.avg_image_time_s[1], old_data.avg_image_time_s[1]],
            color = [colors["JuliaMSI"], colors["JuliaMSIBulk"], colors["julia_mzML_imzML"]],
            width = bar_width
        )
        
        # Plot loading memory
        barplot!(ax3, [i + x_offsets[1], i + x_offsets[2], i + x_offsets[3]],
            [new_data.load_memory_mb[1], bulk_data.load_memory_mb[1], old_data.load_memory_mb[1]],
            color = [colors["JuliaMSI"], colors["JuliaMSIBulk"], colors["julia_mzML_imzML"]],
            width = bar_width
        )
        
        # Plot image generation memory
        barplot!(ax4, [i + x_offsets[1], i + x_offsets[2], i + x_offsets[3]],
            [new_data.avg_image_memory_mb[1], bulk_data.avg_image_memory_mb[1], old_data.avg_image_memory_mb[1]],
            color = [colors["JuliaMSI"], colors["JuliaMSIBulk"], colors["julia_mzML_imzML"]],
            width = bar_width
        )
        
        # Calculate speedup factors (higher is better for new library)
        if nrow(new_data) == 1 && nrow(old_data) == 1
            total_time_speedup_new[i] = old_data.total_time_s[1] / new_data.total_time_s[1]
            image_time_speedup_new[i] = old_data.avg_image_time_s[1] / new_data.avg_image_time_s[1]
            load_mem_speedup_new[i] = old_data.load_memory_mb[1] / new_data.load_memory_mb[1]
            img_mem_speedup_new[i] = old_data.avg_image_memory_mb[1] / new_data.avg_image_memory_mb[1]
        end
        if nrow(bulk_data) == 1 && nrow(old_data) == 1
            total_time_speedup_bulk[i] = old_data.total_time_s[1] / bulk_data.total_time_s[1]
            image_time_speedup_bulk[i] = old_data.avg_image_time_s[1] / bulk_data.avg_image_time_s[1]
            load_mem_speedup_bulk[i] = old_data.load_memory_mb[1] / bulk_data.load_memory_mb[1]
            img_mem_speedup_bulk[i] = old_data.avg_image_memory_mb[1] / bulk_data.avg_image_memory_mb[1]
        end

        # Collect data for the summary table CSV
        current_case_name = unique_names[i]
        case_row = df[(df.name .== current_case_name) .& (df.library .== "JuliaMSI"), :] # Any library row will have file_size_mb
        file_size = nrow(case_row) == 1 ? case_row.file_size_mb[1] : NaN

        push!(summary_df_rows, Dict(
            "DS" => unique_names[i],
            "File Size (MB)" => file_size,
            "JMSI-TT Speedup" => round(total_time_speedup_new[i], digits=2),
            "JMSI-IT Speedup" => round(image_time_speedup_new[i], digits=2),
            "JMSI-LM Speedup" => round(load_mem_speedup_new[i], digits=2),
            "JMSI-IM Speedup" => round(img_mem_speedup_new[i], digits=2),
            "JBulk-TT Speedup" => round(total_time_speedup_bulk[i], digits=2),
            "JBulk-IT Speedup" => round(image_time_speedup_bulk[i], digits=2),
            "JBulk-LM Speedup" => round(load_mem_speedup_bulk[i], digits=2),
            "JBulk-IM Speedup" => round(img_mem_speedup_bulk[i], digits=2)
        ))
    end
    
    # Create DataFrame from collected rows
    summary_df = DataFrame(summary_df_rows)
    
    # Save the summary table to CSV
    CSV.write(joinpath(RESULTS_DIR, "benchmark_summary_table.csv"), summary_df)
    println("\nBenchmark summary table saved to $(joinpath(RESULTS_DIR, "benchmark_summary_table.csv"))")

    # Add overall title
    Label(fig1[1, :], "MSI Library Performance Benchmark", fontsize = 72, font = :bold)
    
    # Add a legend
    Legend(fig1[3, :], # Moved to row 3, spanning all columns
        [PolyElement(color = colors["JuliaMSI"]), PolyElement(color = colors["JuliaMSIBulk"]), PolyElement(color = colors["julia_mzML_imzML"])],
        [library_labels["JuliaMSI"], library_labels["JuliaMSIBulk"], library_labels["julia_mzML_imzML"]],
        "Library",
        orientation = :horizontal,
        tellwidth = false, tellheight = true,
        halign = :center, valign = :top, # Centered horizontally, aligned to top
        margin = (10, 10, 10, 10)
    )
    
    # Adjust layout spacing
    colgap!(g1, 1, 30)
    rowgap!(g1, 1, 20)
    
    # Ensure all rows and columns in g1 expand to fill available space
    rowsize!(g1, 1, Auto())
    rowsize!(g1, 2, Auto())
    colsize!(g1, 1, Auto())
    colsize!(g1, 2, Auto())
    
    # Ensure the row containing the legend expands in the top-level figure layout
    rowsize!(fig1.layout, 3, Auto())

    

    # Save plots

    save(joinpath(RESULTS_DIR, "performance_dashboard.png"), fig1, px_per_unit = 1)    
    println("\nEnhanced visualizations saved to $RESULTS_DIR")
end

function plot_walkthrough_results(walkthrough_df)
    fig = Figure(size=(3600, 2400))
    
    # Plot 1: Total time per file
    ax1 = Axis(fig[1, 1], 
               title="Time per File (Load + Generate 1 Image)",
               ylabel="Time (seconds)",
               xticks=(collect(1:nrow(walkthrough_df)), walkthrough_df.name),
               xticklabelrotation=π/8)
    
    # Grouped bars
    barplot!(ax1, collect(1:nrow(walkthrough_df)) .- 0.2, walkthrough_df.JuliaMSI_time_s, 
             width=0.4, color="#2E86AB", label="JuliaMSI (New)")
    barplot!(ax1, collect(1:nrow(walkthrough_df)) .+ 0.2, walkthrough_df.Old_time_s,
             width=0.4, color="#A23B72", label="julia_mzML_imzML (Old)")
    
    # Plot 2: Speedup factor
    ax2 = Axis(fig[1, 2],
               title="Speedup Factor (Old/New Time)",
               ylabel="Speedup (Higher = Better)",
               xticks=(collect(1:nrow(walkthrough_df)), walkthrough_df.name),
               xticklabelrotation=π/8)
    
    # Bars above 1 = new library is faster
    colors = [x >= 1 ? "#2CA02C" : "#E53935" for x in walkthrough_df.speedup]
    barplot!(ax2, collect(1:nrow(walkthrough_df)), walkthrough_df.speedup, color=colors)
    hlines!(ax2, [1.0], color=:black, linestyle=:dash, linewidth=2)
    
    # Plot 3: Memory comparison
    ax3 = Axis(fig[2, 1],
               title="Memory Usage per File",
               ylabel="Memory (MB)",
               xticks=(collect(1:nrow(walkthrough_df)), walkthrough_df.name),
               xticklabelrotation=π/8)
    
    barplot!(ax3, collect(1:nrow(walkthrough_df)) .- 0.2, walkthrough_df.JuliaMSI_mem_mb,
             width=0.4, color="#2E86AB", label="JuliaMSI")
    barplot!(ax3, collect(1:nrow(walkthrough_df)) .+ 0.2, walkthrough_df.Old_mem_mb,
             width=0.4, color="#A23B72", label="Old Library")
    
    # Plot 4: Time vs File Size (scatter)
    ax4 = Axis(fig[2, 2],
               title="Processing Time vs File Size",
               xlabel="File Size (MB)",
               ylabel="Time (seconds)")
    
    scatter!(ax4, walkthrough_df.file_size_mb, walkthrough_df.JuliaMSI_time_s,
             color="#2E86AB", markersize=20, label="JuliaMSI")
    scatter!(ax4, walkthrough_df.file_size_mb, walkthrough_df.Old_time_s,
             color="#A23B72", markersize=20, label="Old Library")
    
    # Add linear trend lines
    if nrow(walkthrough_df) > 1
        linreg_new = hcat(fill(1.0, nrow(walkthrough_df)), walkthrough_df.file_size_mb) \ walkthrough_df.JuliaMSI_time_s
        linreg_old = hcat(fill(1.0, nrow(walkthrough_df)), walkthrough_df.file_size_mb) \ walkthrough_df.Old_time_s
        
        x_vals = [minimum(walkthrough_df.file_size_mb), maximum(walkthrough_df.file_size_mb)]
        lines!(ax4, x_vals, linreg_new[1] .+ linreg_new[2] .* x_vals, color="#2E86AB", linewidth=2)
        lines!(ax4, x_vals, linreg_old[1] .+ linreg_old[2] .* x_vals, color="#A23B72", linewidth=2)
    end
    
    # Add legends
    Legend(fig[0, :], [PolyElement(color="#2E86AB"), PolyElement(color="#A23B72")],
           ["JuliaMSI (New)", "julia_mzML_imzML (Old)"], 
           orientation=:horizontal, tellwidth=false, framevisible=true)
    
    # Adjust layout
    rowgap!(fig.layout, 1)
    colgap!(fig.layout, 1)
    
    save(joinpath(RESULTS_DIR, "multi_file_walkthrough.png"), fig, px_per_unit=1)
    println("\nMulti-file walkthrough plot saved.")
end

function plot_real_world_pipeline_results(pipeline_df::DataFrame)
    fig = Figure(size=(3600, 2400)) # Larger figure for more detailed plots

    # Filter out rows with Inf or NaN values in critical plotting columns
    # Create a copy to avoid modifying the original DataFrame
    filtered_df = deepcopy(pipeline_df)
    
    # Filter for finite values in key metrics for plotting
    filter!(row -> isfinite(row.JuliaMSI_time_s) && isfinite(row.Old_time_s) && 
                   isfinite(row.speedup) && isfinite(row.JuliaMSI_mem_mb) && 
                   isfinite(row.Old_mem_mb) && isfinite(row.memory_reduction_pct) &&
                   isfinite(row.slices_per_second_new) && isfinite(row.slices_per_second_old), 
            filtered_df)

    if isempty(filtered_df)
        println("No finite data points for real-world pipeline plotting after filtering. Skipping plot generation.")
        save(joinpath(RESULTS_DIR, "real_world_pipeline.png"), fig, px_per_unit=1)
        return
    end

    num_rows_filtered = nrow(filtered_df)
    x_positions = collect(1:num_rows_filtered)
    
    # Plot 1: Total Time Comparison
    ax1 = Axis(fig[1, 1],
               title="Real-World Pipeline: Total Time",
               ylabel="Time (seconds)",
               xticks=(x_positions, filtered_df.name),
               xticklabelrotation=π/8)
    
    barplot!(ax1, x_positions .- 0.2, filtered_df.JuliaMSI_time_s,
             width=0.4, color="#2E86AB", label="JuliaMSI (New)")
    barplot!(ax1, x_positions .+ 0.2, filtered_df.Old_time_s,
             width=0.4, color="#A23B72", label="Old Library")
    
    # Plot 2: Speedup Factor
    ax2 = Axis(fig[1, 2],
               title="Real-World Pipeline: Speedup Factor (Old/New Time)",
               ylabel="Speedup (Higher = Better)",
               xticks=(x_positions, filtered_df.name),
               xticklabelrotation=π/8)
    
    colors_speedup = [x >= 1 ? "#2CA02C" : "#E53935" for x in filtered_df.speedup]
    barplot!(ax2, x_positions, filtered_df.speedup, color=colors_speedup)
    hlines!(ax2, [1.0], color=:black, linestyle=:dash, linewidth=2)
    
    # Plot 3: Memory Usage Comparison
    ax3 = Axis(fig[2, 1],
               title="Real-World Pipeline: Memory Usage",
               ylabel="Memory (MB)",
               xticks=(x_positions, filtered_df.name),
               xticklabelrotation=π/8)
    
    barplot!(ax3, x_positions .- 0.2, filtered_df.JuliaMSI_mem_mb,
             width=0.4, color="#2E86AB", label="JuliaMSI")
    barplot!(ax3, x_positions .+ 0.2, filtered_df.Old_mem_mb,
             width=0.4, color="#A23B72", label="Old Library")

    # Plot 4: Memory Reduction Percentage
    ax4 = Axis(fig[2, 2],
               title="Real-World Pipeline: Memory Reduction vs Old (%)",
               ylabel="Memory Reduction (%)",
               xticks=(x_positions, filtered_df.name),
               xticklabelrotation=π/8)
    
    colors_mem_red = [x > 0 ? "#2CA02C" : "#E53935" for x in filtered_df.memory_reduction_pct]
    barplot!(ax4, x_positions, filtered_df.memory_reduction_pct, color=colors_mem_red)
    hlines!(ax4, [0.0], color=:black, linestyle=:dash, linewidth=2)

    # Plot 5: Slices Per Second
    ax5 = Axis(fig[3, 1],
               title="Real-World Pipeline: Slices per Second",
               ylabel="Slices/sec (Higher = Better)",
               xticks=(x_positions, filtered_df.name),
               xticklabelrotation=π/8)
    
    barplot!(ax5, x_positions .- 0.2, filtered_df.slices_per_second_new,
             width=0.4, color="#2E86AB", label="JuliaMSI")
    barplot!(ax5, x_positions .+ 0.2, filtered_df.slices_per_second_old,
             width=0.4, color="#A23B72", label="Old Library")

    # Overall Legend
    Legend(fig[0, :], [PolyElement(color="#2E86AB"), PolyElement(color="#A23B72")],
           ["JuliaMSI (New)", "julia_mzML_imzML (Old)"],
           orientation=:horizontal, tellwidth=false, framevisible=true)
    
    # Adjust layout
    rowgap!(fig.layout, 1)
    colgap!(fig.layout, 1)
    
    save(joinpath(RESULTS_DIR, "real_world_pipeline.png"), fig, px_per_unit=1)
    println("\nReal-world pipeline plot saved.")
end

# --- Execute ---
mkpath(RESULTS_DIR)
run_all_benchmarks(true) # true: skip processing files # false: create new matrix
