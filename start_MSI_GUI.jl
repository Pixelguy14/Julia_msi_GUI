using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.update()
# Pkg.add(url="https://github.com/CINVESTAV-LABI/julia_mzML_imzML") # with this we ensure it uses the latest library iteration
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
    end
end
