using Pkg
sTime = time()

# Always activate the repository's own project, not the caller's cwd.
Pkg.activate(@__DIR__)

# Keep the user's environment choice. Normal GUI launches use production mode:
# dev mode starts Revise/GenieAutoReload watchers, which adds startup cost and
# duplicate "Watching"/"Autoreloading" log lines. Use JULIA_MSI_DEV=1 to opt in.
if !haskey(ENV, "GENIE_ENV")
    ENV["GENIE_ENV"] = get(ENV, "JULIA_MSI_DEV", "0") == "1" ? "dev" : "prod"
end

# NativeFileDialog uses GTK on Linux. A C/POSIX locale produces the noisy
# "Locale not supported by C library" warning when the dialog opens.
if Sys.islinux() && get(ENV, "LC_ALL", "") == "" &&
   get(ENV, "LANG", "") in ("", "C", "POSIX")
    ENV["LC_ALL"] = "C.UTF-8"
end

# `instantiate` is idempotent and is required even when Manifest.toml exists:
# a manifest can be present while packages are missing from the depot.
try
    Pkg.instantiate()
catch e
    @warn "Pkg.instantiate() failed; retrying after Pkg.resolve()" exception=e
    Pkg.resolve()
    Pkg.instantiate()
end

# Fail before Genie loads if generated assets cannot be written.
include(joinpath(@__DIR__, "src", "Paths.jl"))
ensure_public_dir()

using Logging

# Genie/EzXML prints one warning per Vue/Quasar custom tag while parsing the
# Stipple layout. Those warnings are harmless but make the console unreadable.
# Keep real warnings/errors; allow `JULIA_MSI_VERBOSE=1` to see everything.
struct JuliaMSIFilterLogger <: Logging.AbstractLogger
    inner::Logging.AbstractLogger
    verbose::Bool
end

Logging.min_enabled_level(l::JuliaMSIFilterLogger) =
    Logging.min_enabled_level(l.inner)

Logging.shouldlog(l::JuliaMSIFilterLogger, level, _module, group, id) =
    Logging.shouldlog(l.inner, level, _module, group, id)

Logging.catch_exceptions(l::JuliaMSIFilterLogger) =
    Logging.catch_exceptions(l.inner)

function Logging.handle_message(l::JuliaMSIFilterLogger, level, message, _module,
                                group, id, file, line; kwargs...)
    if !l.verbose
        text = try
            string(message)
        catch
            ""
        end
        # EzXML/libxml2 logs one "Tag <x> invalid" warning per Vue/Quasar custom
        # tag while Genie parses the Stipple layout. Genie 5.35 already silences
        # most of these; keep this guard for older/edge render paths.
        if occursin("XMLError:", text) && occursin("HTML parser", text)
            return nothing
        end
        # Safety net for the Stipple watcher `isempty(nothing)` bug (the root
        # cause is fixed in app.jl by emitting an empty sesstoken meta tag).
        if nameof(_module) === :Stipple &&
           occursin("MethodError(iterate, (nothing", text)
            return nothing
        end
    end
    return Logging.handle_message(l.inner, level, message, _module, group, id,
                                  file, line; kwargs...)
end

Logging.global_logger(JuliaMSIFilterLogger(
    Logging.ConsoleLogger(stderr, Logging.Info),
    get(ENV, "JULIA_MSI_VERBOSE", "0") == "1"
))

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

            # Do not remove sessions from another currently-running instance.
            # Only clean directories that have been idle for at least an hour.
            (time() - mtime(full_path)) < 3600 && continue

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


# Load and configure Genie. Do not redirect stderr: errors must remain visible
# to the user and to the log/diagnostics workflow.
Genie.loadapp()

# Start the Genie server
@async begin
    up(host="127.0.0.1", port=1481)
    eTime=round(time()-sTime,digits=3)
    println("Julia MSI GUI took $(eTime) seconds booting")
    url = "http://127.0.0.1:1481"
    # Open the URL in the default web browser based on the OS.
    # `start` and `Start-Process` are shell built-ins, not executables, so
    # Windows must go through `cmd /c`.
    if Sys.isapple()
        @async try
            run(`open $url`)
        catch e
            @warn "Could not open the browser automatically" url=url exception=e
        end
    elseif Sys.islinux()
        @async try
            run(`xdg-open $url`)
        catch e
            @warn "Could not open the browser automatically" url=url exception=e
        end
    elseif Sys.iswindows()
        @async try
            run(`cmd /c start "" $url`)
        catch e
            @warn "Could not open the browser automatically" url=url exception=e
        end
    end
end
wait()
