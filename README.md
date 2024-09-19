# GraphsMatching

[![Build status](https://github.com/JuliaGraphs/GraphsMatching.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaGraphs/GraphsMatching.jl/actions/workflows/ci.yml?query=branch%3Amaster)
[![Code coverage](http://codecov.io/github/JuliaGraphs/GraphsMatching.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaGraphs/GraphsMatching.jl?branch=master)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Matching algorithms on top of [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl).

## Usage

The results of any matching is returned as a `MatchingResult` struct
containing the `mate` and `weight` fields.

### Perfect matching

```julia
g = complete_graph(4)
w = Dict{Edge,Float64}()
w[Edge(1,3)] = 10
w[Edge(1,4)] = 0.5
w[Edge(2,3)] = 11
w[Edge(2,4)] = 2
w[Edge(1,2)] = 100

# find the perfect matching of minimum weight
match = minimum_weight_perfect_matching(g, w, 50)
# match.mate[1] == 4
# match.mate[4] == 1
# match.mate[2] == 3
# match.mate[3] == 2
# match.weight ≈ 11.5
```

### Maximum weight matching

A maximum weight matching is solved as a Linear Programming
problem and requires an LP optimizer for bipartite graphs and a MILP solver for general graphs respecting the [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl) optimizer interface. A list of solvers can be found in the [JuMP documentation](http://www.juliaopt.org/JuMP.jl/v0.19.0/installation/#Getting-Solvers-1).

```julia
using JuMP, Cbc #import a MILP solver
g = complete_graph(3)
w = zeros(3,3)
w[1,2] = 1
w[3,2] = 1
w[1,3] = 1
match = maximum_weight_matching(g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),w)
# match.weight ≈ 1
```
