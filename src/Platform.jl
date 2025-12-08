# src/Platform.jl

module Platform

using Base.Threads # For Threads.nthreads()

# Export functions for system capability detection
export detect_system_capabilities

"""
    detect_system_capabilities()

Detects and returns key system capabilities such as total memory and number of CPU threads.

# Returns
- A `NamedTuple` with fields:
    - `total_memory_gb::Float64`: Total physical memory in gigabytes.
    - `num_cpu_threads::Int`: Number of available CPU threads.
"""
function detect_system_capabilities()::NamedTuple
    total_memory_bytes = Sys.total_memory()
    total_memory_gb = total_memory_bytes / (1024^3) # Convert bytes to gigabytes

    num_cpu_threads = Sys.CPU_THREADS # Total logical CPU cores

    return (total_memory_gb=total_memory_gb, num_cpu_threads=num_cpu_threads)
end

end # module Platform
