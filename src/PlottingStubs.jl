# src/PlottingStubs.jl
#
# These names used to live in MSI_src and were exported from the core module.
# They are GUI/plotting helpers that require CairoMakie, so their
# implementations now live in julia_imzML_visual.jl. Keeping explicit stubs
# here gives downstream code an actionable error instead of an UndefVarError,
# while keeping `using MSI_src` lightweight for headless use.
#
# The stubs are intentionally NOT exported.

const _PLOTTING_MOVED_MESSAGE = """
This function is a GUI/plotting helper and is no longer implemented in the
MSI_src core. It moved to the GUI layer (julia_imzML_visual.jl) so that
`using MSI_src` does not load CairoMakie.

Load the GUI module that defines it, or call the corresponding helper from your
application/plotting layer. See the JuliaMSI README for the migration note.
"""

function plot_slice(args...; kwargs...)
    error(_PLOTTING_MOVED_MESSAGE)
end

function plot_slices(args...; kwargs...)
    error(_PLOTTING_MOVED_MESSAGE)
end

function generate_colorbar_image(args...; kwargs...)
    error(_PLOTTING_MOVED_MESSAGE)
end
