# scripts/precompile.jl
using Pkg
Pkg.activate(dirname(@__DIR__))

using Genie
using HTTP

@info "Precompiling application..."

# Load the app. Mmap=false is recommended for compilation
Genie.loadapp(pwd(); Mmap=false)

@info "Starting server for precompilation."
# Start the server in the background
server = up(1481, "127.0.0.1", async=true)

# Give the server a moment to start up
sleep(20)

base_url = "http://127.0.0.1:1481"
@info "Hitting routes to precompile..."

try
    # Make a GET request to the home page
    response = HTTP.get(base_url * "/")
    @info "GET / -> Status: $(response.status)"

    # === IMPORTANT ===
    # For best performance, you should add more HTTP requests here
    # to hit ALL your important routes and API endpoints.
    # This ensures the code for every page gets compiled.
    # Example:
    # HTTP.get(base_url * "/contact")
    # HTTP.post(base_url * "/submit-data", [], "{\"key\":\"value\"}")

catch e
    @error "Could not hit routes during precompilation." exception=(e, catch_backtrace())

finally
    # Stop the server
    @info "Shutting down server."
    down()
end

@info "Precompilation tracing finished."
