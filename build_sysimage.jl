# build_sysimage.jl
using PackageCompiler

println("--- Starting Custom System Image Build ---")
println("This process can take 5-10 minutes.")

# Define the image path with OS-specific extension
ext = Sys.iswindows() ? ".dll" : Sys.isapple() ? ".dylib" : ".so"
sysimage_path = "MSI_sysimage" * ext

# List dependencies to include in the sysimage for maximum stability
# We have removed graphical heavy-lifters (CairoMakie, PlotlyBase) to ensure 
# the build completes successfully on systems with standard RAM.
packages_to_include = [
    :Genie,
    :GenieFramework,
    :PlotlyBase,
    :Libz,
    :DataFrames,
    :Serialization,
    :Printf,
    :JSON
]

create_sysimage(
    packages_to_include,
    sysimage_path = sysimage_path,
    precompile_execution_file = "precompile_script.jl",
    incremental = true,
    filter_stdlibs = false
)

println("--- System Image Build Successful! ---")
println("To start Julia with this image, use:")
println("julia --threads auto --project=. --sysimage $(sysimage_path) start_MSI_GUI.jl")
