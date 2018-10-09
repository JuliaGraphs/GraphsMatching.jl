# LightGraphsMatching

[![Build Status](https://travis-ci.org/JuliaGraphs/LightGraphsMatching.jl.svg?branch=master)](https://travis-ci.org/JuliaGraphs/LightGraphsMatching.jl)

[![Coverage Status](https://coveralls.io/repos/github/JuliaGraphs/LightGraphsMatching.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGraphs/LightGraphsMatching.jl?branch=master)

[![codecov](https://codecov.io/gh/JuliaGraphs/LightGraphsMatching.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGraphs/LightGraphsMatching.jl)

Matching algorithms on top of [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl).

## Usage

The results of any matching is returned as a `MatchingResult` struct
containing the `mate` and `weight` fields.

### Perfect matching

```julia
g =CompleteGraph(4)
w =Dict{Edge,Float64}()
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
problem and requires a LP solver respecting the [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) solver
interface. See MathProgBase 
[documentation](http://mathprogbasejl.readthedocs.io/en/latest/solvers.html) for more details.

```julia
using Cbc: CbcSolver #import a LP solver
g = CompleteGraph(3)
w = zeros(3,3)
w[1,2] = 1
w[3,2] = 1
w[1,3] = 1
match = maximum_weight_matching(g,CbcSolver(),w)
# match.weight ≈ 1
```
