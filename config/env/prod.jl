# config/env/prod.jl
# This file is loaded when GENIE_ENV is set to "prod"

# Define PLUGINS_WITH_ASSETS to prevent UndefVarError in production mode
# It's typically an empty vector if no specific plugins with assets are tracked.
const PLUGINS_WITH_ASSETS = []
