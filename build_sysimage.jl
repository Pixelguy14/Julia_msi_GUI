# build_sysimage.jl
using PackageCompiler

headless = "--headless" in ARGS

println("--- Starting Custom System Image Build ---")
println("Mode: ", headless ? "Headless (Scripting/Docker)" : "Full GUI (Genie App)")
println("This process can take 5-15 minutes.")

ext = Sys.iswindows() ? ".dll" : Sys.isapple() ? ".dylib" : ".so"
sysimage_name = headless ? "sys_msi_headless" * ext : "sys_msi_gui" * ext
sysimage_path = joinpath(@__DIR__, sysimage_name)

packages_to_include = Symbol[
    :MSI_src, # <--- Bake our core library
    :DataFrames,
    :SparseArrays,
    :Mmap,
    :Libz,
    :Serialization,
    :Printf,
    :JSON
]

if !headless
    append!(packages_to_include, [
        :Genie,
        :GenieFramework,
        :PlotlyBase
    ])
end

create_sysimage(
    packages_to_include,
    sysimage_path = sysimage_path,
    precompile_execution_file = "precompile_script.jl",
    incremental = true,
    filter_stdlibs = false
)

println("--- System Image Build Successful! ---")
println("Created: $(sysimage_path)")
if headless
    println("To use in a headles script:")
    println("julia --threads auto --project=. --sysimage $(sysimage_name) my_script.jl")
else
    println("To start Julia with the GUI image, use:")
    println("julia --threads auto --project=. --sysimage $(sysimage_name) start_MSI_GUI.jl")
end
