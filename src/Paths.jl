# src/Paths.jl
# Shared filesystem locations for JuliaMSI.
#
# All application-owned files live under the project's `public/` directory.
# Keeping those paths in one place avoids mixing `pwd()`-relative and
# `@__DIR__`-relative paths (which is especially fragile on Windows).

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const PUBLIC_DIR = joinpath(PROJECT_ROOT, "public")
const DEFAULT_REGISTRY_PATH = joinpath(PUBLIC_DIR, "registry.json")

"""
    ensure_public_dir(; path::AbstractString = PUBLIC_DIR) -> String

Ensures that the application's `public/` directory exists and is writable.

JuliaMSI writes generated bitmaps, masks, and `registry.json` into this
directory, so a read-only or missing directory is a fatal configuration error.
A short write probe is used instead of relying only on permission bits, because
Windows ACLs and read-only attributes are not reliably reflected by `iswritable`.
"""
function ensure_public_dir(; path::AbstractString = PUBLIC_DIR)
    try
        mkpath(path)
    catch e
        error("JuliaMSI cannot create the public directory at \"$path\". " *
              "Move the project to a writable location or fix the directory permissions. " *
              "Original error: $(sprint(showerror, e))")
    end

    isdir(path) ||
        error("JuliaMSI expected a directory at \"$path\", but it is not a directory.")

    probe = joinpath(path, ".juliamsi_write_probe_$(getpid())")
    try
        open(probe, "w") do io
            write(io, "ok")
        end
    catch e
        error("JuliaMSI cannot write to \"$path\". The application requires write " *
              "access to its public directory for generated images, masks, and " *
              "registry.json. Move the project to a writable location or fix the " *
              "directory permissions. Original error: $(sprint(showerror, e))")
    finally
        isfile(probe) && rm(probe; force=true)
    end

    return String(path)
end
