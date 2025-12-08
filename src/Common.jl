# src/Common.jl - Updated with BloomFilter
using Base.Threads

# --- Buffer Pooling ---
"""
    BufferPool

A thread-safe pool of `Vector{UInt8}` buffers, categorized by their size.
This helps reduce memory allocations and garbage collection overhead by reusing buffers.

# Fields
- `lock`: A `ReentrantLock` to ensure thread-safe access to the pool.
- `buffers`: A dictionary mapping buffer size (in bytes) to a list of available buffers of that size.
"""
mutable struct BufferPool
    lock::ReentrantLock
    buffers::Dict{Int, Vector{Vector{UInt8}}}

    function BufferPool()
        new(ReentrantLock(), Dict{Int, Vector{Vector{UInt8}}}())
    end
end

"""
    get_buffer!(pool::BufferPool, size::Int) -> Vector{UInt8}

Retrieves a `Vector{UInt8}` buffer of at least `size` bytes from the pool.
If no suitable buffer is available, a new one is allocated.

# Arguments
- `pool`: The `BufferPool` to retrieve the buffer from.
- `size`: The minimum desired size of the buffer in bytes.

# Returns
- A `Vector{UInt8}` buffer.
"""
function get_buffer!(pool::BufferPool, size::Int)::Vector{UInt8}
    lock(pool.lock) do
        # Find an existing buffer that is large enough
        for (buffer_size, buffer_list) in pool.buffers
            if buffer_size >= size && !isempty(buffer_list)
                buffer = pop!(buffer_list)
                # Resize if necessary (and if it's significantly larger than needed)
                if length(buffer) > size * 2 # Heuristic: if buffer is more than twice the requested size
                    return Vector{UInt8}(undef, size) # Allocate new smaller buffer
                else
                    return buffer
                end
            end
        end
        # No suitable buffer found, allocate a new one
        return Vector{UInt8}(undef, size)
    end
end

"""
    release_buffer!(pool::BufferPool, buffer::Vector{UInt8})

Returns a `Vector{UInt8}` buffer to the pool for future reuse.

# Arguments
- `pool`: The `BufferPool` to return the buffer to.
- `buffer`: The `Vector{UInt8}` buffer to release.
"""
function release_buffer!(pool::BufferPool, buffer::Vector{UInt8})
    lock(pool.lock) do
        buffer_size = length(buffer)
        if !haskey(pool.buffers, buffer_size)
            pool.buffers[buffer_size] = Vector{Vector{UInt8}}()
        end
        push!(pool.buffers[buffer_size], buffer)
    end
    return nothing
end


"""
    SimpleBufferPool

A dramatically simplified buffer pool that avoids complex locking.
"""
mutable struct SimpleBufferPool
    buffers::Dict{Int, Vector{Vector{UInt8}}}
    max_pool_size::Int
end

SimpleBufferPool() = SimpleBufferPool(Dict{Int, Vector{Vector{UInt8}}}(), 50)

function get_buffer!(pool::SimpleBufferPool, size::Int)::Vector{UInt8}
    # Check for existing buffers of exact size first
    if haskey(pool.buffers, size) && !isempty(pool.buffers[size])
        return pop!(pool.buffers[size])
    end
    
    # No suitable buffer found, allocate new one
    return Vector{UInt8}(undef, size)
end

function release_buffer!(pool::SimpleBufferPool, buffer::Vector{UInt8})
    size = length(buffer)
    
    if !haskey(pool.buffers, size)
        pool.buffers[size] = Vector{Vector{UInt8}}()
    end
    
    # Limit pool size to prevent memory bloat
    if length(pool.buffers[size]) < pool.max_pool_size
        push!(pool.buffers[size], buffer)
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
const DEFAULT_CACHE_SIZE = 100
const DEFAULT_NUM_BINS = 2000
const DEFAULT_BLOOM_FILTER_SIZE = 10000  # Default size for empty spectra
const DEFAULT_FALSE_POSITIVE_RATE = 0.01

"""
    validate_spectrum_data(mz, intensity, id)

Checks a spectrum for common data integrity issues.
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
"""
function empty_bloom_filter()::BloomFilter{Int}
    return BloomFilter{Int}(size=DEFAULT_BLOOM_FILTER_SIZE, hash_count=3)
end

"""
    AtomicFlag

A thread-safe boolean flag using atomic operations.
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
