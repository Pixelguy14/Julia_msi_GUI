# Respect an explicitly selected environment. This makes
# `GENIE_ENV=prod julia --project=. start_MSI_GUI.jl` actually work;
# previously this file unconditionally forced dev mode.
if !haskey(ENV, "GENIE_ENV")
    ENV["GENIE_ENV"] = "prod"
end
