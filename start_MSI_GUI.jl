using Pkg
sTime=time()
Pkg.activate(".")
Pkg.instantiate()
Pkg.gc()

using Genie

# Load and configure Genie
Genie.loadapp()

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