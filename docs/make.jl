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
