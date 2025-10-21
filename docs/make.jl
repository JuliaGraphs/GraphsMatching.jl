using Documenter
using GraphsMatching
import Graphs

makedocs(
    modules = [GraphsMatching],
    format = Documenter.HTML(),
    sitename = "GraphsMatching",
	pages = Any[
		"Getting started"         => "index.md",
		"API"                     => "API.md",
		"Internals"               => "internals.md",
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
