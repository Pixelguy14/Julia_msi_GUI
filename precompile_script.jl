# precompile_script.jl
# This script exercises the core MSI processing kernels (imzML and mzML) 
# to ensure they are precompiled into the custom system image.

using MSI_src
using Mmap
using Base64
using PlotlyBase

# --- DYNAMIC WARM-UP WORKLOAD --- #

function create_warmup_datasets(dir, T_mz::Type, T_int::Type)
    # 1. Create minimal imzML/ibd
    imzml_path = joinpath(dir, "warmup_$(T_mz)_$(T_int).imzML")
    ibd_path = joinpath(dir, "warmup_$(T_mz)_$(T_int).ibd")
    
    n_points = 10
    n_pixels = 1
    mz_bytes = n_points * sizeof(T_mz)
    int_bytes = n_points * sizeof(T_int)
    
    # Write some dummy binary data (m/z must be sorted for validation)
    open(ibd_path, "w") do f
        write(f, sort(rand(T_mz, n_points)))
        write(f, rand(T_int, n_points))
    end
    
    # Minimal but VALID imzML structure that passes axes_config_img
    imzml_content = """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzML" version="1.1.0">
  <cvList count="2">
    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="3.30.0" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/ms/mzML/ontology/psi-ms.obo"/>
    <cv id="IMS" fullName="Imaging MS Ontology" version="0.9.1" URI="http://www.imzml.org/ontology/imzML1.1.0.obo"/>
  </cvList>
  <fileDescription>
    <fileContent>
      <cvParam cvRef="IMS" accession="IMS:1000031" name="processed" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="universally unique identifier" value="00000000-0000-0000-0000-000000000000"/>
      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd MD5 HTTP" value="00000000000000000000000000000000"/>
    </fileContent>
  </fileDescription>
  <referenceableParamGroupList count="2">
    <referenceableParamGroup id="mz_array_settings">
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
      <cvParam cvRef="MS" accession="MS:$(T_mz == Float64 ? "1000523" : "1000521")" name="$(T_mz == Float64 ? "64-bit float" : "32-bit float")" value=""/>
      <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensity_array_settings">
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
      <cvParam cvRef="MS" accession="MS:$(T_int == Float64 ? "1000523" : "1000521")" name="$(T_int == Float64 ? "64-bit float" : "32-bit float")" value=""/>
      <cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" value=""/>
    </referenceableParamGroup>
  </referenceableParamGroupList>
  <run id="run_1">
    <spectrumList count="1">
      <spectrum index="0" id="scan=1" defaultArrayLength="$n_points">
        <cvParam cvRef="IMS" accession="IMS:1000050" name="coordinate x" value="1"/>
        <cvParam cvRef="IMS" accession="IMS:1000051" name="coordinate y" value="1"/>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mz_array_settings"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external data" value=""/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external offset" value="0"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external array length" value="$n_points"/>
            <cvParam cvRef="IMS" accession="IMS:1000101" name="external encoded length" value="$mz_bytes"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="intensity_array_settings"/>
            <cvParam cvRef="IMS" accession="IMS:1000102" name="external data" value=""/>
            <cvParam cvRef="IMS" accession="IMS:1000103" name="external offset" value="$mz_bytes"/>
            <cvParam cvRef="IMS" accession="IMS:1000104" name="external array length" value="$n_points"/>
            <cvParam cvRef="IMS" accession="IMS:1000101" name="external encoded length" value="$int_bytes"/>
            <binary/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>
"""
    write(imzml_path, imzml_content)

    # 2. Create minimal mzML (Base64)
    mzml_path = joinpath(dir, "warmup_$(T_mz)_$(T_int).mzML")
    b64_mz = base64encode(sort(rand(T_mz, n_points)))
    b64_int = base64encode(rand(T_int, n_points))
    
    mzml_content = """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzML" version="1.1.0">
  <run id="run_1">
    <spectrumList count="1">
      <spectrum index="0" id="scan=1" defaultArrayLength="$n_points">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="$(length(b64_mz))">
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
            <cvParam cvRef="MS" accession="MS:$(T_mz == Float64 ? "1000523" : "1000521")" name="" value=""/>
            <binary>$b64_mz</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="$(length(b64_int))">
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
            <cvParam cvRef="MS" accession="MS:$(T_int == Float64 ? "1000523" : "1000521")" name="" value=""/>
            <binary>$b64_int</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
  <indexList count="1"><index name="spectrum"><offset idRef="scan=1">PLACEHOLDER</offset></index></indexList>
  <indexListOffset>ID_OFFSET</indexListOffset>
</mzML>
"""
    spec_start = findfirst("<spectrum", mzml_content).start - 1
    index_start = findfirst("<indexList", mzml_content).start - 1
    mzml_content = replace(mzml_content, "PLACEHOLDER" => string(spec_start))
    mzml_content = replace(mzml_content, "ID_OFFSET" => string(index_start))
    
    write(mzml_path, mzml_content)
    
    return imzml_path, mzml_path
end

println("Starting robust precompilation warmup (imzML + mzML)...")

try
    mktempdir() do tmp_dir
        for (T_mz, T_int) in [(Float32, Float32), (Float64, Float32)]
            println("Exercising pathways for $T_mz and $T_int...")
            imzml_path, mzml_path = create_warmup_datasets(tmp_dir, T_mz, T_int)
            
            try
                # Exercise imzML pathway
                imzml_data = MSI_src.load_imzml_lazy(imzml_path; use_mmap=true)
                MSI_src.get_multiple_mz_slices(imzml_data, [10.0, 20.0], 0.1)
                
                # Exercise Sprint 2 Streaming Pipeline
                config = MSI_src.PipelineConfig(
                    steps=[
                        MSI_src.StreamingStep(:baseline_correction, Dict(:method => :snip, :iterations => 5)),
                        MSI_src.StreamingStep(:normalization, Dict(:method => :tic)),
                        MSI_src.StreamingStep(:peak_picking, Dict(:method => :profile, :snr_threshold => 3.0))
                    ],
                    num_bins=2000
                )
                try
                    MSI_src.process_dataset!(imzml_data, config)
                catch
                    # Ignore matrix mismatch errors from tiny random data 
                end
                
                # Exercise mzML pathway
                mzml_data = MSI_src.load_mzml_lazy(mzml_path)
                MSI_src.GetSpectrum(mzml_data, 1)
            catch e
                @warn "Warmup failed for $T_mz/$T_int: $e"
            end
        end

        # --- COMMON KERNELS ---
        println("Precompiling image processing kernels...")
        dummy_slice = rand(Float32, 4, 4)
        MSI_src.TrIQ(dummy_slice, 256, 1.0)
        MSI_src.save_bitmap(joinpath(tmp_dir, "d.bmp"), zeros(UInt8, 4, 4), MSI_src.ViridisPalette)
        
        # Plotly Warmup
        p = PlotlyBase.Plot(PlotlyBase.scatter(x=1:2, y=[1,2]))
    end
catch e
    if !(e isa InterruptException)
        @warn "Global Warmup failed: $e"
    end
end
println("Precompilation workload completed.")
