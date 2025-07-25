"""
    minimum_weight_perfect_matching(g, w::Dict{Edge,Real}; tmaxscale)
    minimum_weight_perfect_matching(g, w::Dict{Edge,Real}, cutoff; tmaxscale)
    minimum_weight_perfect_matching(g, w::Dict{Edge,Real}, algorithm::AbstractMinimumWeightPerfectMatchingAlgorithm; tmaxscale)
    minimum_weight_perfect_matching(g, w::Dict{Edge,Real}, cutoff, algorithm::AbstractMinimumWeightPerfectMatchingAlgorithm; tmaxscale)

Given a graph `g` and an edgemap `w` containing weights associated to edges,
returns a matching with the mimimum total weight among the ones containing
exactly `nv(g)/2` edges.

Edges in `g` not present in `w` will not be considered for the matching.

You can use the `algorithm` argument to specify the algorithm to use.

A `cutoff` argument can be given, to reduce the computational time by
excluding edges with weights higher than the cutoff (effective only for some algorithms,
not for the default `LEMONMWPMAlgorithm`).

When the weights are non-integer types, the keyword argument `tmaxscale` can be used to
scale the weights to integer values.
In case of error try to change `tmaxscale` (default is `tmaxscale=10`).
The scaling is as follows:
```
tmax = typemax(Int32) / tmaxscale
weight = round(Int32, (weight-minimum_weight) / max(maximum_weight-minimum_weight, 1) * tmax)
```

The returned object is of type [`MatchingResult`](@ref).

See also: [`BlossomVAlgorithm`](@ref)
"""
function minimum_weight_perfect_matching end

"""
    AbstractMinimumWeightPerfectMatchingAlgorithm

Abstract type that allows users to pass in their preferred algorithm

See also: [`minimum_weight_perfect_matching`](@ref), [`BlossomVAlgorithm`](@ref)
"""
abstract type AbstractMinimumWeightPerfectMatchingAlgorithm end

"""
    BlossomVAlgorithm()

Use the BlossomV algorithm to find the minimum weight perfect matching.
Depends on the BlossomV.jl package. You have to call `using BlossomV`
before using this algorithm.

This algorithm dispatches to the BlossomV library, a C library for finding
minimum weight perfect matchings in general graphs.
The BlossomV library is not open source, and thus we can not distribute a
a precompiled binary with GraphsMatching.jl. We attempt to build it on
installation, but unlike typical Julia packages, the build process is prone
to failure if you do not have all the dependencies installed.
If BlossomV.jl does not work on your system,
consider using the LEMONGraphs.jl algorithm instead (the default algorithm),
which we distribute precompiled on all platforms.

See also: [`minimum_weight_perfect_matching`](@ref), [`LEMONMWPMAlgorithm`](@ref)
"""
struct BlossomVAlgorithm <: AbstractMinimumWeightPerfectMatchingAlgorithm end

"""
    LEMONMWPMAlgorithm()

Use the LEMON C++ implementation of minimum weight perfect matching.

See also: [`minimum_weight_perfect_matching`](@ref), [`BlossomVAlgorithm`](@ref)
"""
struct LEMONMWPMAlgorithm <: AbstractMinimumWeightPerfectMatchingAlgorithm end

function minimum_weight_perfect_matching(
    g::Graph, w::Dict{E,U}
) where {U<:Integer,E<:Edge}
    return minimum_weight_perfect_matching(g, w, LEMONMWPMAlgorithm())
end

function minimum_weight_perfect_matching(
    g::Graph, w::Dict{E,U}, cutoff::Real, algorithm::AbstractMinimumWeightPerfectMatchingAlgorithm=LEMONMWPMAlgorithm(); kws...
) where {U<:Real,E<:Edge}
    wnew = Dict{E,U}()
    for (e, c) in w
        if c <= cutoff
            wnew[e] = c
        end
    end
    return minimum_weight_perfect_matching(g, wnew, algorithm; kws...)
end

function minimum_weight_perfect_matching(
    g::Graph, w::Dict{E,U}, algorithm::AbstractMinimumWeightPerfectMatchingAlgorithm=LEMONMWPMAlgorithm(); tmaxscale=10.0
) where {U<:AbstractFloat,E<:Edge}
    wnew = Dict{E,Int32}()
    cmax = maximum(values(w))
    cmin = minimum(values(w))

    tmax = typemax(Int32) / tmaxscale # /10 is kinda arbitrary,
    # hopefully high enough to not occur in overflow problems
    for (e, c) in w
        wnew[e] = round(Int32, (c - cmin) / max(cmax - cmin, 1) * tmax)
    end
    match = minimum_weight_perfect_matching(g, wnew, algorithm)
    weight = zero(U)
    for i in 1:nv(g)
        j = match.mate[i]
        if j > i
            weight += get(w, E(i, j), zero(U))
        end
    end
    return MatchingResult(weight, match.mate)
end

function minimum_weight_perfect_matching(g::Graph, w::Dict{E,U}, ::LEMONMWPMAlgorithm) where {U<:Integer,E<:Edge}
    max = 2*abs(maximum(values(w)))
    weights = [-get(w, e, max) for e in edges(g)]
    totweight, mate = LEMONGraphs.maxweightedperfectmatching(g, weights)
    return MatchingResult(-totweight, mate)
end
