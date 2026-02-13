
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MSI_src
using Test
using Base.Threads

# --- Test Configuration ---
const TEST_FILE = "/home/pixel/Documents/Cinvestav_2025/Analisis/Stomach/Stomach_DHB_uncompressed.imzML"

println("Starting concurrency test with file: $TEST_FILE")

if !isfile(TEST_FILE)
    error("Test file not found: $TEST_FILE")
end

# --- Simulation of the Concurrency Issue ---

# Global "UI State" variable simulates the app's msi_data
global_msi_data = nothing
const GLOBAL_LOCK = ReentrantLock()

function ui_load_file(path)
    global global_msi_data
    lock(GLOBAL_LOCK) do
        if global_msi_data !== nothing
            close(global_msi_data)
        end
        println("UI: Loading file...")
        global_msi_data = OpenMSIData(path)
        println("UI: File loaded.")
    end
end

function ui_close_file()
    global global_msi_data
    lock(GLOBAL_LOCK) do
        if global_msi_data !== nothing
            println("UI: Closing file...")
            close(global_msi_data)
            global_msi_data = nothing
            println("UI: File closed.")
        end
    end
end

# Mock pipeline function that mimics `run_full_pipeline` in app.jl
# Crucially, it uses its own local `pipeline_msi_data` as per the fix
function run_mock_pipeline_isolated(path)
    println("Pipeline: Starting isolated pipeline...")
    
    # 1. Open ISOLATED instance
    pipeline_msi_data = OpenMSIData(path)
    println("Pipeline: Opened isolated MSIData instance.")
    
    try
        # 2. Simulate reading spectra
        indices = 1:min(100, length(pipeline_msi_data.spectra_metadata)) # Read first 100 spectra
        
        # Artificial delay to allow "user" interaction
        sleep(0.5) 
        
        println("Pipeline: Reading spectra...")
        Threads.@threads for i in indices
            # Use local instance
            mz, int = GetSpectrum(pipeline_msi_data, i)
            # Simulate processing work
            sum(int) 
        end
        println("Pipeline: Finished reading spectra successfully.")
        return true
    catch e
        println("Pipeline: CRASHED with error: $e")
        return false
    finally
        close(pipeline_msi_data)
        println("Pipeline: Closed isolated MSIData instance.")
    end
end

# Mock pipeline that uses GLOBAL instance (The BUGGY version)
function run_mock_pipeline_buggy()
    println("Buggy Pipeline: Starting...")
    # Uses global_msi_data directly
    
    try
        global global_msi_data
        if global_msi_data === nothing
            println("Buggy Pipeline: No data loaded!")
            return false
        end
        
        local_ref = global_msi_data # Still points to same object
        
        indices = 1:min(100, length(local_ref.spectra_metadata))
        
        sleep(0.5)
        
        println("Buggy Pipeline: Reading spectra from shared object...")
        Threads.@threads for i in indices
            # This will fail if ui_close_file() happens concurrently
            mz, int = GetSpectrum(local_ref, i)
            sum(int)
        end
        println("Buggy Pipeline: Success (Unexpected if concurrency worked)")
        return true
    catch e
        println("Buggy Pipeline: CRASHED as expected: $e")
        return false
    end
end


# --- execute checks ---

@testset "Concurrency Crash Fix Verification" begin
    
    # 1. Setup: Load file initially
    ui_load_file(TEST_FILE)
    
    # 2. Test the FIX: Isolated Pipeline
    println("\n--- Testing Fixed (Isolated) Pipeline ---")
    
    t_pipeline = Threads.@spawn run_mock_pipeline_isolated(TEST_FILE)
    
    # Simulate user closing/reloading file while pipeline runs
    sleep(0.2) 
    ui_close_file()
    
    # Wait for pipeline
    success = fetch(t_pipeline)
    @test success == true
    println("Fixed pipeline result: ", success ? "PASSED" : "FAILED")


    # 3. Test the BUG: Global Pipeline (Optional, to prove it crashes without fix)
    # Uncomment to verify the bug exists if needed, but we assume it does based on user report.
    # println("\n--- Testing Buggy (Shared) Pipeline ---")
    # ui_load_file(TEST_FILE)
    # t_buggy = Threads.@spawn run_mock_pipeline_buggy()
    # sleep(0.2)
    # ui_close_file()
    # buggy_success = fetch(t_buggy)
    # println("Buggy pipeline result: ", buggy_success ? "PASSED (No crash?)" : "FAILED (Crashed as expected)")
    
end
