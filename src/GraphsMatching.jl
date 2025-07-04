module GraphsMatching

using Graphs

using SparseArrays: spzeros

using JuMP
using MathOptInterface
const MOI = MathOptInterface
import LEMONGraphs
using Hungarian

export MatchingResult,
    maximum_weight_matching,
    maximum_weight_matching_reduction,
    maximum_weight_maximal_matching,
    minimum_weight_perfect_matching,
    HungarianAlgorithm,
    LPAlgorithm,
    BlossomVAlgorithm,
    LEMONMWPMAlgorithm

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if BlossomVAlgorithm in argtypes
                print(io,"""\nPlease first import `BlossomV` to make BlossomV-based matching available.""")
            end
        end
    end
end

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

include("lp.jl")
include("maximum_weight_matching.jl")
include("minimum_weight_perfect_matching.jl")
include("hungarian.jl")
include("maximum_weight_maximal_matching.jl")

end # module
