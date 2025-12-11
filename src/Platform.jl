# src/Platform.jl

using Base.Threads # For Threads.nthreads()

"""
    PlatformProfile

A struct containing a snapshot of the system's platform characteristics.

# Fields
- `device_type::Symbol`: Categorization of the device (e.g., `:server`, `:desktop`, `:laptop`, `:mobile`).
- `os::Symbol`: Operating system (e.g., `:linux`, `:windows`, `:macos`).
- `num_cpu_threads::Int`: Number of logical CPU cores.
- `total_ram_gb::Float64`: Total physical memory in gigabytes.
"""
struct PlatformProfile
    device_type::Symbol
    os::Symbol
    num_cpu_threads::Int
    total_ram_gb::Float64
end

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

"""
    detect_platform_profile()::PlatformProfile

Detects and returns a `PlatformProfile` for the current system.

# Returns
- A `PlatformProfile` struct containing system characteristics.
"""
function detect_platform_profile()::PlatformProfile
    system_caps = detect_system_capabilities()
    total_ram_gb = system_caps.total_memory_gb
    num_cpu_threads = system_caps.num_cpu_threads

    # Determine OS
    os_symbol = if Sys.islinux()
        :linux
    elseif Sys.iswindows()
        :windows
    elseif Sys.isapple()
        :macos
    else
        :unknown
    end

    # Heuristic for device type
    device_type = :unknown
    if total_ram_gb >= 16 && num_cpu_threads >= 8
        device_type = :desktop
        if total_ram_gb >= 64 && num_cpu_threads >= 32
            device_type = :server
        end
    elseif total_ram_gb < 16 || num_cpu_threads < 4
        device_type = :laptop # Or mobile if even smaller, but focusing on typical Julia dev machines
    end

    return PlatformProfile(device_type, os_symbol, num_cpu_threads, total_ram_gb)
end
