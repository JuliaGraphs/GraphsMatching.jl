module GraphsMatching

using Graphs

using SparseArrays: spzeros

using JuMP
using MathOptInterface
const MOI = MathOptInterface
using BlossomV: BlossomV # 'using BlossomV'  leads to naming conflicts with JuMP
using Hungarian

export MatchingResult,
    weight,
    is_matched_vertex,
    matching_vertex,
    matching_vertices,
    matched_edges,
    maximum_weight_matching,
    maximum_weight_maximal_matching,
    minimum_weight_perfect_matching,
    HungarianAlgorithm,
    LPAlgorithm

"""
    struct MatchingResult{U}
        weight::U
        mate::Vector{Int}
    end

A type representing the result of a matching algorithm.

    weight: total weight of the matching

    mate:    `mate[i] = j` if vertex `i` is matched to vertex `j`.
             `mate[i] = -1` for unmatched vertices.
"""
struct MatchingResult{U<:Real}
    weight::U
    mate::Vector{Int}
end

"""
    weight(m::MatchingResult)

    returns the total weight of the matching
"""
weight(m::MatchingResult) = m.weight

"""
    is_matched_vertex(m::MatchingResult, u)

    returns true if vertex `u` is matched with another vertex
"""
is_matched_vertex(m::MatchingResult, u) = (m.mate[u] != -1)

"""
    matching_vertex(m::MatchingResult, u)

    returns the vertex matched with `u` (-1 if `u` is not matched)
"""
matching_vertex(m::MatchingResult, u) = m.mate[u]

"""
    matching_vertices(m::MatchingResult)

    returns a list of the matching vertices (-1 if the vertex is not matched)
"""
matching_vertices(m::MatchingResult) = m.mate

"""
    matching_vertex(m::MatchingResult, u)

    returns a list of the matched edges
"""
matched_edges(m::MatchingResult) =
    [Edge(u, v) for (u, v) in enumerate(matching_vertices(m)) if (u < v)]

include("lp.jl")
include("maximum_weight_matching.jl")
include("blossomv.jl")
include("hungarian.jl")
include("maximum_weight_maximal_matching.jl")

end # module
