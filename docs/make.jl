using Documenter
using LightGraphsMatching
import LightGraphs; const lg = LightGraphs

makedocs(
    modules = [LightGraphsMatching],
    format = :html,
    sitename = "LightGraphsMatching",
	pages = Any[
		"Getting started"         => "index.md",
    ]
)

deploydocs(
    deps        = nothing,
    make        = nothing,
    repo        = "github.com/JuliaGraphs/LightGraphsMatching.jl.git",
    target      = "build",
    julia       = "0.6",
    osname      = "linux"
)
