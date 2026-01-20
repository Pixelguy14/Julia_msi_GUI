using Pkg
sTime=time()
Pkg.activate(".")
ENV["GENIE_ENV"] = "dev"
#ENV["GENIE_ENV"] = "prod"

manifest_path = joinpath(@__DIR__, "Manifest.toml")

# Only instantiate in development mode
if get(ENV, "GENIE_ENV", "dev") != "prod" || !isfile(manifest_path)
    @info "Development environment detected. Instantiating packages..."
    
    if !isfile(manifest_path)
        @info "Manifest.toml not found. Generating it based on Project.toml..."
    end

    Pkg.resolve()
    Pkg.instantiate()
    Pkg.gc()
end

using Genie

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