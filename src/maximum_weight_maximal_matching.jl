"""
    AbstractMaximumWeightMaximalMatchingAlgorithm
Abstract type that allows users to pass in their preferred algorithm
"""
abstract type AbstractMaximumWeightMaximalMatchingAlgorithm end

"""
    LPAlgorithm <: AbstractMaximumWeightMaximalMatchingAlgorithm
Forces the maximum_weight_maximal_matching function to use a linear programming formulation.
"""
struct LPAlgorithm <: AbstractMaximumWeightMaximalMatchingAlgorithm end

function maximum_weight_maximal_matching(
    g::Graph, 
    w::AbstractMatrix{T}, 
    algorithm::LPAlgorithm, 
    solver = nothing
) where {T<:Real}
    if ! isa(solver, AbstractMathProgSolver)
        error("The keyword argument solver must be an AbstractMathProgSolver, as accepted by JuMP.")
    end

    return maximum_weight_maximal_matching_lp(g, solver, w)
end

"""
    HungarianAlgorithm <: AbstractMaximumWeightMaximalMatchingAlgorithm
Forces the maximum_weight_maximal_matching function to use the Hungarian algorithm.
"""
struct HungarianAlgorithm <: AbstractMaximumWeightMaximalMatchingAlgorithm end

function maximum_weight_maximal_matching(
    g::Graph, 
    w::AbstractMatrix{T}, 
    algorithm::HungarianAlgorithm, 
    solver = nothing
) where {T<:Real}
    return maximum_weight_maximal_matching_hungarian(g, w)
end

"""
    maximum_weight_maximal_matching{T<:Real}(g, w::Dict{Edge,T})

Given a bipartite graph `g` and an edge map `w` containing weights associated to edges,
returns a matching with the maximum total weight among the ones containing the
greatest number of edges.

Edges in `g` not present in `w` will not be considered for the matching.

A `cutoff` keyword argument can be given, to reduce computational times
excluding edges with weights lower than the cutoff.

Finally, a specific algorithm can be chosen (`algorithm` keyword argument); 
each algorithm has specific dependencies. For instance: 

- If `algorithm=HungarianAlgorithm()` (the default), the package Hungarian.jl is used. 
  This algorithm is always polynomial in time, with complexity O(n³). 
- If `algorithm=LPAlgorithm()`, the package JuMP.jl and one of its supported solvers is required. 
  In this case, the algorithm relies on a linear relaxation on of the matching problem, which is
  guaranteed to have integer solution on bipartite graphs. A solver must be provided with 
  the `solver` keyword parameter. 

The returned object is of type `MatchingResult`.
"""
function maximum_weight_maximal_matching(
        g::Graph, 
        w::AbstractMatrix{T}; 
        cutoff = nothing,
        algorithm::AbstractMaximumWeightMaximalMatchingAlgorithm = HungarianAlgorithm(), 
        solver = nothing
    ) where {T<:Real}

    if cutoff != nothing && ! isa(cutoff, Real)
        error("The cutoff value must be of type Real or nothing.")
    end

    if cutoff != nothing
        return maximum_weight_maximal_matching(g, cutoff_weights(w, cutoff), algorithm, solver)
    else
        return maximum_weight_maximal_matching(g, w, algorithm, solver)
    end
end

"""
    cutoff_weights copies the weight matrix with all elements below cutoff set to 0
"""
function cutoff_weights(w::AbstractMatrix{T}, cutoff::R) where {T<:Real, R<:Real}
    wnew = copy(w)
    for j in 1:size(w,2)
        for i in 1:size(w,1)
            if wnew[i,j] < cutoff
                wnew[i,j] = zero(T)
            end
        end
    end
    wnew
end

@deprecate maximum_weight_maximal_matching(g::Graph, solver::AbstractMathProgSolver, w::AbstractMatrix{T}, cutoff::R) where {T<:Real, R<:Real} maximum_weight_maximal_matching(g, w, algorithm=LPAlgorithm(), cutoff=cutoff, solver=solver)
