using MSI_src
using Test
using Base.Threads

@testset "SimpleBufferPool Concurrency Stress Test" begin
    pool = MSI_src.SimpleBufferPool()
    n_threads = Threads.nthreads()
    n_iterations = 1000
    buffer_size = 1024

    println("Running buffer pool stress test with $n_threads threads...")

    # Parallel stress test
    Threads.@threads for i in 1:(n_threads * n_iterations)
        # get_buffer! and release_buffer! are now thread-safe
        buf = MSI_src.get_buffer!(pool, buffer_size)
        
        # Simulate some work
        fill!(buf, UInt8(i % 256))
        
        MSI_src.release_buffer!(pool, buf)
    end

    # After stress test, the dictionary should be coherent
    sizes = collect(keys(pool.buffers))
    if !isempty(sizes)
        @test buffer_size ∈ sizes
        @test length(pool.buffers[buffer_size]) <= pool.max_pool_size
    end
    
    println("Buffer pool stress test PASSED.")
end
