using Test
using MSI_src

@testset "JuliaMSI Core Systems" begin
    @testset "Basic Instantiation" begin
        # 1. Test basic metadata structure initialization
        meta = SpectrumMetadata(
            Int32(1), Int32(2), 
            "test_id", 
            :sample, 
            CENTROID, 
            SpectrumAsset(Float64, false, Int64(0), 100, :mz, 0.0, 0.0), 
            SpectrumAsset(Float32, true, Int64(100), 50, :intensity, 0.0, 0.0)
        )
        @test meta.x == 1
        @test meta.y == 2
        @test meta.mode == CENTROID
        @test meta.mz_asset.format == Float64
        @test meta.int_asset.format == Float32
        @test meta.int_asset.is_compressed == true

        # 2. Test cache pool and MSIData structural stability
        source = MzMLSource([], Float64, Float32, nothing) # Empty source for structural testing
        
        # Test constructor doesn't throw
        # Constructor signature: (source, metadata, instrument_meta, dims, coordinate_map, cache_size)
        msi_data = MSIData(
            source, 
            [meta], 
            nothing, # instrument_meta
            (100, 100), # dims
            nothing, # coord map
            10 # cache size
        )
        
        @test msi_data.image_dims == (100, 100)
        @test length(msi_data.spectra_metadata) == 1
        @test msi_data.cache_size == 10
    end

    @testset "Buffer & Cache Subsystems" begin
        pool = SimpleBufferPool()
        # Test basic retrieval
        buf = get_buffer!(pool, 1024)
        @test length(buf) == 1024
        
        # Test release
        release_buffer!(pool, buf)
        @test length(pool.buffers[1024]) == 1
        
        # Test reuse
        buf2 = get_buffer!(pool, 1024)
        @test buf === buf2 # Should return the EXACT same buffer object
    end

    @testset "Bloom Filter Hashing and Probing" begin
        bf = MSI_src.BloomFilter{Int}(1000, 0.01)
        inserted = collect(1:1000)
        foreach(x -> push!(bf, x), inserted)

        # Bloom filters must never produce false negatives.
        @test all(x -> x in bf, inserted)
        @test all(x -> MSI_src.contains(bf, x), inserted)
        @test MSI_src.false_positive_rate(bf) > 0.0
        @test MSI_src.fill_ratio(bf) > 0.0
    end

    @testset "Application Paths" begin
        # `PUBLIC_DIR` must be absolute and independent of the caller's cwd.
        @test isabspath(MSI_src.PUBLIC_DIR)
        @test MSI_src.ensure_public_dir() == MSI_src.PUBLIC_DIR
    end
end
