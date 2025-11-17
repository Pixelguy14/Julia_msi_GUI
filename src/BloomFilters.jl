# src/BloomFilter.jl

"""
    BloomFilter{T}

A simple, efficient Bloom filter implementation for probabilistic set membership testing.

# Fields
- `bits::BitVector`: The underlying bit array
- `size::Int`: Number of bits in the filter
- `hash_count::Int`: Number of hash functions to use
- `seed::UInt64`: Random seed for hash functions
- `count::Int`: Number of elements added (for monitoring)
"""
mutable struct BloomFilter{T}
    bits::BitVector
    size::Int
    hash_count::Int
    seed::UInt64
    count::Int
end

"""
    StreamingBloomFilter

A memory-efficient Bloom filter that doesn't store all values at once.
"""
mutable struct StreamingBloomFilter
    bits::BitVector
    size::Int
    hash_count::Int
    seed::UInt64
    count::Int
    # Streaming state
    current_chunk::Vector{Int}
    chunk_size::Int
end

"""
    BloomFilter{T}(expected_elements::Int, false_positive_rate::Float64=0.01; kwargs...) -> Return type

Creates a Bloom filter optimized for the expected number of elements and desired false positive rate.

# Arguments

- `expected_elements::Int`: Argument description
- `false_positive_rate::Float64`: Argument description
    (**Default**: `0.01`)

# Keywords

- `seed::Union{UInt32,UInt64}`: Keyword description
    (**Default**: `0x12345678`)
"""
function BloomFilter{T}(expected_elements::Int, false_positive_rate::Float64=0.01; seed::Union{UInt32,UInt64}=0x12345678) where T
    # Convert seed to UInt64 for consistency
    seed_uint64 = UInt64(seed)
    
    # Calculate optimal parameters using standard formulas
    size = optimal_bit_size(expected_elements, false_positive_rate)
    hash_count = optimal_hash_count(expected_elements, size)
    
    bits = falses(size)
    return BloomFilter{T}(bits, size, hash_count, seed_uint64, 0)
end

"""
    optimal_bit_size(n::Int, p::Float64) -> Int

Calculates the optimal number of bits for a Bloom filter

# Arguments

- `n::Int`: Expected number of elements
- `p::Float64`: Desired false positive rate
"""
function optimal_bit_size(n::Int, p::Float64)::Int
    if p <= 0.0 || p >= 1.0
        throw(ArgumentError("False positive rate must be between 0 and 1"))
    end
    # m = - (n * ln(p)) / (ln(2)^2)
    m = ceil(Int, - (n * log(p)) / (log(2)^2))
    return max(m, 1)
end

"""
    optimal_hash_count(n::Int, m::Int) -> Int

Calculates the optimal number of hash functions for a Bloom filter.
    
# Arguments

- `n::Int`: Expected number of elements
- `p::Float64`: Desired false positive rate
"""
function optimal_hash_count(n::Int, m::Int)::Int
    if n <= 0 || m <= 0
        return 1
    end
    # k = (m / n) * ln(2)
    k = max(1, round(Int, (m / n) * log(2)))
    return min(k, 8)  # Practical limit to avoid too many hashes
end

"""
    hash_functions(item::T, count::Int, size::Int, seed::UInt64) -> Vector{Int}

Generates multiple hash values for an item using double hashing technique.
"""
function hash_functions(item::T, count::Int, size::Int, seed::UInt64) where T
    # Use Julia's built-in hash with different seeds
    hashes = Vector{Int}(undef, count)
    
    # First hash with the main seed
    h1 = hash(item, seed)
    h2 = hash(item, seed + 1)
    
    for i in 1:count
        # Double hashing: h_i = h1 + i * h2
        combined_hash = UInt64(h1) + UInt64(i) * UInt64(h2)
        hashes[i] = (combined_hash % UInt64(size)) + 1  # 1-based indexing
    end
    
    return hashes
end

"""
    Base.push!(bf::BloomFilter{T}, item::T)

Adds an element to the Bloom filter.

# Arguments

- `bf::BloomFilter{T}`: Argument description
- `item::T`: Argument description
"""
function Base.push!(bf::BloomFilter{T}, item::T) where T
    hashes = hash_functions(item, bf.hash_count, bf.size, bf.seed)
    
    for h in hashes
        bf.bits[h] = true
    end
    bf.count += 1
    
    return bf
end

"""
    Base.in(item::T, bf::BloomFilter{T}) -> Bool

Checks whether an element is possibly in the Bloom filter.
Returns `true` if the element might be in the set, `false` if it's definitely not.
"""
function Base.in(item::T, bf::BloomFilter{T})::Bool where T
    hashes = hash_functions(item, bf.hash_count, bf.size, bf.seed)
    
    for h in hashes
        if !bf.bits[h]
            return false  # Definitely not in the set
        end
    end
    return true  # Possibly in the set
end

"""
    contains(bf::BloomFilter{T}, item::T) -> Bool

Alias for `in()` for compatibility with your existing code.
"""
contains(bf::BloomFilter{T}, item::T) where T = item in bf

"""
    add!(bf::BloomFilter{T}, item::T)

Alias for `push!()` for compatibility with your existing code.
"""
add!(bf::BloomFilter{T}, item::T) where T = push!(bf, item)

"""
    false_positive_rate(bf::BloomFilter) -> Float64

Estimates the current false positive rate of the Bloom filter.
"""
function false_positive_rate(bf::BloomFilter)::Float64
    if bf.count == 0
        return 0.0
    end
    # Theoretical false positive rate: (1 - e^(-k * n / m)) ^ k
    k = bf.hash_count
    n = bf.count
    m = bf.size
    return (1 - exp(-k * n / m)) ^ k
end

"""
    fill_ratio(bf::BloomFilter)::Float64 return count(bf.bits) / length(bf.bits) end -> Return type

Returns the fraction of bits that are set to 1.

# Arguments

- `bf::BloomFilter`: Argument description
"""
function fill_ratio(bf::BloomFilter)::Float64
    return count(bf.bits) / length(bf.bits)
end

"""
    is_empty(bf::BloomFilter)::Bool return bf.count == 0 end -> Return type

Checks if the Bloom filter is empty (no elements added).

# Arguments

- `bf::BloomFilter`: Argument description
"""
function is_empty(bf::BloomFilter)::Bool
    return bf.count == 0
end

"""
    reset!(bf::BloomFilter) fill!(bf.bits, false) bf.count = 0 return bf end -> Return type

Clears the Bloom filter, removing all elements.

# Arguments

- `bf::BloomFilter`: Argument description
"""
function reset!(bf::BloomFilter)
    fill!(bf.bits, false)
    bf.count = 0
    return bf
end

# Specialized constructor for empty Bloom filters
"""
    BloomFilter{T}(; kwargs...) -> Return type

Description of the function

# Keywords

- `size::Int`: Keyword description
    (**Default**: `100`)
- `hash_count::Int`: Keyword description
    (**Default**: `3`)
- `seed::Union{UInt32,UInt64}`: Keyword description
    (**Default**: `0x12345678`)
"""
function BloomFilter{T}(;size::Int=100, hash_count::Int=3, seed::Union{UInt32,UInt64}=0x12345678) where T
    seed_uint64 = UInt64(seed)
    bits = falses(size)
    return BloomFilter{T}(bits, size, hash_count, seed_uint64, 0)
end
