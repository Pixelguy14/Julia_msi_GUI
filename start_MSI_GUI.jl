using Pkg
sTime=time()
Pkg.activate(".")
ENV["GENIE_ENV"] = "dev"
#ENV["GENIE_ENV"] = "prod"

manifest_path = joinpath(@__DIR__, "Manifest.toml")

# Selective instantiation for faster startup
if get(ENV, "GENIE_ENV", "dev") != "prod" && !isfile(manifest_path)
    @info "Development environment detected and Manifest.toml missing. Instantiating packages..."
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.gc()
elseif get(ENV, "GENIE_ENV", "dev") != "prod"
    @info "Manifest.toml found. Skipping Pkg.instantiate() for faster boot. Delete Manifest.toml if you need to re-instantiate."
    ENV["GENIE_ENV"] = "prod"
end

using Genie

# --- Cross-Platform Startup Cleanup ---
# Remove orphaned GenieSessionFileSession directories from previous runs.
# These accumulate in the OS temp directory as jl_XXXXXX folders containing
# serialized session files (64-char hex filenames). Over long sessions or
# after crashes, they can consume gigabytes of disk space.
function cleanup_orphaned_sessions()
    tmp = Base.tempdir()
    cleaned_count = 0
    cleaned_bytes = 0

    for entry in readdir(tmp; join=false)
        # Only target directories matching Julia's temp naming pattern
        startswith(entry, "jl_") || continue
        full_path = joinpath(tmp, entry)
        isdir(full_path) || continue

        # Validate: a Genie session dir contains files with 64-char hex names
        try
            contents = readdir(full_path)
            isempty(contents) && continue

            # Check if at least one file matches the 64-char hex session ID pattern
            is_session_dir = any(contents) do f
                length(f) == 64 && all(c -> c in "0123456789abcdef", f)
            end
            is_session_dir || continue

            # Safe to remove — this is an orphaned Genie session directory
            dir_size = sum(filesize(joinpath(full_path, f)) for f in contents; init=0)
            rm(full_path; recursive=true, force=true)
            cleaned_count += 1
            cleaned_bytes += dir_size
        catch e
            @debug "Skipping $entry during cleanup: $e"
        end
    end

    if cleaned_count > 0
        size_mb = round(cleaned_bytes / (1024^2), digits=1)
        @info "Startup cleanup: removed $cleaned_count orphaned session dir(s), freed $(size_mb) MB"
    end
end

cleanup_orphaned_sessions()


# Load and configure Genie
Genie.loadapp()

# Remove html parser error discrepancy
redirect_stderr(devnull)

# Start the Genie server
@async begin
    up(host="127.0.0.1", port=1481)
    eTime=round(time()-sTime,digits=3)
    println("Julia MSI GUI took $(eTime) seconds booting")
    url = "http://127.0.0.1:1481"
    # Open the URL in the default web browser based on the OS
    if Sys.isapple()
        @async run(`open $url`) # For macOS
    elseif Sys.islinux()
        @async run(`xdg-open $url`) # For Linux
    elseif Sys.iswindows()
        @async run(`start $url`) # For Windows
        @async run(`explorer $url`)
        @async run(`Start-Process $url`) # For Windows
    end
end
wait()