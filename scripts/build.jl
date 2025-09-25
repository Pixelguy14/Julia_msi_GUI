# scripts/build.jl
using PackageCompiler
using Pkg

project_dir = dirname(@__DIR__)
Pkg.activate(project_dir)

# Get all the direct dependencies from Project.toml
deps = keys(Pkg.project().dependencies)
sysimage_path = joinpath(project_dir, "build", "JuliaMSI_sysimage.so")
precompile_script = joinpath(project_dir, "scripts", "precompile.jl")

@info "Starting system image compilation. This may take several minutes..."

create_sysimage(
    deps;
    sysimage_path=sysimage_path,
    precompile_execution_file=precompile_script
)

@info "Compilation finished!"
@info "You can now run the fast-starting application using the following command:"
@info "julia --project=. -J build/JuliaMSI_sysimage.so start_MSI_GUI.jl"
