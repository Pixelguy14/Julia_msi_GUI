# src/ResourcePool.jl
using Base.Threads

"""
    ResourcePool{T}

A thread-safe pool for reusing objects of type `T` to minimize allocations and GC pressure.
Specifically designed for high-performance computing tasks where large buffers are needed
repeatedly across multiple threads.

# Fields:
- `pool::Vector{T}`: The underlying storage for idle resources.
- `lock::ReentrantLock`: Ensures thread-safe access to the pool.
- `max_size::Int`: Maximum number of resources to hold in the pool.
- `constructor::Function`: A function to create a new resource if the pool is empty.
"""
mutable struct ResourcePool{T}
    pool::Vector{T}
    lock::ReentrantLock
    max_size::Int
    constructor::Function
end

"""
    aligned_vector(::Type{T}, n::Int; alignment::Int=64) where T

Creates a `Vector{T}` that is aligned to `alignment` bytes.
Note: In modern Julia, standard vectors are often 16 or 64 byte aligned, but for
HPC we ensure this by allocating slightly more and using a view, or using
specific pointers. For simplicity and performance, we use a small hack:
allocating a larger array and taking a 64-byte aligned view.
"""
function aligned_vector(::Type{T}, n::Int; alignment::Int=64) where T
    # Allocate enough space to find an aligned starting point
    raw = Vector{UInt8}(undef, n * sizeof(T) + alignment)
    ptr = Int(pointer(raw))
    off = (alignment - (ptr % alignment)) % alignment
    # Return a reinterpret view of the aligned segment
    return reinterpret(T, view(raw, (off + 1):(off + n * sizeof(T))))
end

"""
    ResourcePool{T}(constructor::Function; max_size::Int=2 * nthreads())

Creates a new `ResourcePool` for resources of type `T`.
"""
function ResourcePool{T}(constructor::Function; max_size::Int=2 * nthreads()) where T
    return ResourcePool{T}(T[], ReentrantLock(), max_size, constructor)
end

"""
    acquire(pool::ResourcePool{T}) -> T

Retrieves a resource from the pool. If the pool is empty, a new resource is created
using the constructor.
"""
function acquire(pool::ResourcePool{T}) where T
    lock(pool.lock) do
        if !isempty(pool.pool)
            return pop!(pool.pool)
        end
    end
    # Create new resource outside of lock to minimize contention
    return pool.constructor()
end

"""
    release!(pool::ResourcePool{T}, resource::T)

Returns a resource to the pool for later reuse. If the pool is already at `max_size`,
the resource is allowed to be garbage collected.
"""
function release!(pool::ResourcePool{T}, resource::T) where T
    lock(pool.lock) do
        if length(pool.pool) < pool.max_size
            push!(pool.pool, resource)
        end
    end
    return nothing
end

"""
    with_resource(f::Function, pool::ResourcePool{T})

Acquires a resource from the pool, executes the function `f(resource)`, and
automatically releases the resource back to the pool when finished.
"""
function with_resource(f::Function, pool::ResourcePool{T}) where T
    resource = acquire(pool)
    try
        return f(resource)
    finally
        release!(pool, resource)
    end
end
