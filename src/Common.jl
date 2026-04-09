using Base.Threads
using Mmap

# POSIX madvise constants
const MADV_NORMAL = 0
const MADV_RANDOM = 1
const MADV_SEQUENTIAL = 2
const MADV_WILLNEED = 3
const MADV_DONTNEED = 4

"""
    posix_madvise(buffer::AbstractArray, advice::Integer)

A safe wrapper for the OS `madvise` system call. Signals the kernel about the 
access pattern for a memory-mapped region. Currently supports Linux/Unix systems.
"""
function posix_madvise(buffer::AbstractArray, advice::Integer)
    @static if Sys.isunix()
        try
            ptr = pointer(buffer)
            len = sizeof(buffer)
            # ccall(:madvise, return_type, (arg_types...), args...)
            ret = ccall(:madvise, Int32, (Ptr{Cvoid}, Csize_t, Int32), ptr, len, Int32(advice))
            return ret == 0
        catch
            return false
        end
    else
        return false # Not supported on this OS
    end
end

# --- Buffer Pooling ---



"""
    SimpleBufferPool

A dramatically simplified buffer pool that avoids complex locking.

# Arguments:
- `buffers::Dict{Int, Vector{Vector{UInt8}}}`: A dictionary mapping buffer size (in bytes) to a list of available buffers of that size.
- `max_pool_size::Int`: The maximum number of buffers to keep in the pool.
- `lock::ReentrantLock`: A lock to ensure thread-safe access to the pool.
"""
mutable struct SimpleBufferPool
    buffers::Dict{Int, Vector{Vector{UInt8}}}
    max_pool_size::Int
    lock::ReentrantLock  # Add lock for thread safety
end

"""
    SimpleBufferPool()

Creates a new `SimpleBufferPool` with a default maximum pool size of 50.
"""
SimpleBufferPool() = SimpleBufferPool(Dict{Int, Vector{Vector{UInt8}}}(), 50, ReentrantLock())

"""
    get_buffer!(pool::SimpleBufferPool, size::Int) -> Vector{UInt8}

Retrieves a `Vector{UInt8}` buffer of at least `size` bytes from the pool.
If no suitable buffer is available, a new one is allocated.

# Arguments:
- `pool`: The `SimpleBufferPool` to retrieve the buffer from.
- `size`: The minimum desired size of the buffer in bytes.

# Returns:
- A `Vector{UInt8}` buffer.
"""
function get_buffer!(pool::SimpleBufferPool, size::Int)::Vector{UInt8}
    buf = lock(pool.lock) do
        # Check for existing buffers of exact size first
        if haskey(pool.buffers, size) && !isempty(pool.buffers[size])
            return pop!(pool.buffers[size])
        end
        return nothing
    end
    
    if buf !== nothing
        return buf
    end
    
    # No suitable buffer found, allocate new one (outside lock to reduce contention)
    return Vector{UInt8}(undef, size)
end

"""
    release_buffer!(pool::SimpleBufferPool, buffer::Vector{UInt8})

Returns a `Vector{UInt8}` buffer to the pool for future reuse.

# Arguments:
- `pool`: The `SimpleBufferPool` to return the buffer to.
- `buffer`: The `Vector{UInt8}` buffer to release.
"""
function release_buffer!(pool::SimpleBufferPool, buffer::Vector{UInt8})
    size = length(buffer)
    
    lock(pool.lock) do
        if !haskey(pool.buffers, size)
            pool.buffers[size] = Vector{Vector{UInt8}}()
        end
        
        # Limit pool size to prevent memory bloat
        if length(pool.buffers[size]) < pool.max_pool_size
            push!(pool.buffers[size], buffer)
        end
    end
    # If pool is full, let buffer get GC'd
end

# --- Unified Error Types ---
abstract type MSIError <: Exception end

struct FileFormatError <: MSIError
    msg::String
end

struct SpectrumNotFoundError <: MSIError
    id
end

struct InvalidSpectrumError <: MSIError
    id
    reason::String
end

# --- Shared Constants ---
const DEFAULT_CACHE_SIZE = 100 # Default cache size for spectra
const DEFAULT_NUM_BINS = 2000 # Default number of bins for spectra
const DEFAULT_BLOOM_FILTER_SIZE = 10000  # Default size for empty spectra
const DEFAULT_FALSE_POSITIVE_RATE = 0.01 # Default false positive rate for bloom filters

"""
    validate_spectrum_data(mz, intensity, id)

Checks a spectrum for common data integrity issues.

# Arguments:
- `mz`: The m/z values of the spectrum.
- `intensity`: The intensity values of the spectrum.
- `id`: The ID of the spectrum.

# Returns:
- `true` if the spectrum is valid.

# Throws:
- `InvalidSpectrumError` if the spectrum is invalid.
"""
function validate_spectrum_data(mz::AbstractVector, intensity::AbstractVector, id)
    if length(mz) != length(intensity)
        throw(InvalidSpectrumError(id, "m/z and intensity arrays have different lengths"))
    end

    if !all(isfinite, mz) || !all(isfinite, intensity)
        throw(InvalidSpectrumError(id, "Spectrum contains NaN or Inf values"))
    end

    if !issorted(mz)
        throw(InvalidSpectrumError(id, "m/z array is not sorted"))
    end

    return true
end

"""
    create_bloom_filter_for_spectrum(mz::AbstractVector{<:Real}) -> BloomFilter{Int}

Memory-optimized version that processes data in chunks to reduce temporary allocations
and avoid holding large intermediate arrays.

# Arguments:
- `mz`: The m/z values of the spectrum.

# Returns:
- A `BloomFilter{Int}` for the spectrum.
"""
function create_bloom_filter_for_spectrum(mz::AbstractVector{<:Real})::BloomFilter{Int}
    if isempty(mz)
        return empty_bloom_filter()
    end
    
    # Size Bloom filter appropriately
    n_points = length(mz)
    bf = BloomFilter{Int}(n_points, DEFAULT_FALSE_POSITIVE_RATE)
    
    discretization_factor = 100.0
    chunk_size = 10000  # Process in chunks to reduce memory pressure
    
    # Process spectrum data in chunks to avoid large temporary arrays
    for chunk_start in 1:chunk_size:n_points
        chunk_end = min(chunk_start + chunk_size - 1, n_points)
        
        @inbounds for i in chunk_start:chunk_end
            discretized_val = round(Int, mz[i] * discretization_factor)
            push!(bf, discretized_val)
        end
    end
    
    return bf
end

"""
    empty_bloom_filter() -> BloomFilter{Float64}

Creates an empty Bloom filter for spectra with no data.

# Returns:
- An empty `BloomFilter{Float64}`.
"""
function empty_bloom_filter()::BloomFilter{Int}
    return BloomFilter{Int}(size=DEFAULT_BLOOM_FILTER_SIZE, hash_count=3)
end

"""
    AtomicFlag

A thread-safe boolean flag using atomic operations.

# Fields:
- `value::Base.Threads.Atomic{Int}`: The atomic value of the flag.

# Arguments:
- `value::Base.Threads.Atomic{Int}`: The atomic value of the flag.
"""
struct AtomicFlag
    value::Base.Threads.Atomic{Int}
end

AtomicFlag() = AtomicFlag(Base.Threads.Atomic{Int}(0))

function set!(flag::AtomicFlag)
    Base.Threads.atomic_add!(flag.value, 1)
end

function is_set(flag::AtomicFlag)::Bool
    return Base.Threads.atomic_add!(flag.value, 0) > 0
end

function reset!(flag::AtomicFlag)
    Base.Threads.atomic_xchg!(flag.value, 0)
end
