using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.gc()

"""
packages = [
    "GenieFramework", "Libz", "PlotlyBase", "CairoMakie", "Colors", 
    "Statistics", "NaturalSort",  "Genie", 
    "Images", "LinearAlgebra", "NativeFileDialog", "StipplePlotly"
]

# Check for missing packages
for pkg in packages
    if !(pkg in keys(Pkg.dependencies()))
        Pkg.add(pkg)
    end
end

# Add library for mzML imzML from GitHub if missing
if !("julia_mzML_imzML" in keys(Pkg.dependencies()))
    Pkg.add(url="https://github.com/CINVESTAV-LABI/julia_mzML_imzML")
end
"""

using Genie

# Load and configure Genie
Genie.loadapp()

# Start the Genie server
@async begin
    up(host="127.0.0.1", port=1481)
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