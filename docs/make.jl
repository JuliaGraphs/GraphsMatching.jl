using Documenter
using GraphsMatching
import Graphs

makedocs(
    modules = [GraphsMatching],
    format = :html,
    sitename = "GraphsMatching",
	pages = Any[
		"Getting started"         => "index.md",
    ]
)

deploydocs(
    deps        = nothing,
    make        = nothing,
    repo        = "github.com/JuliaGraphs/GraphsMatching.jl.git",
    target      = "build",
    julia       = "0.6",
    osname      = "linux"
)
