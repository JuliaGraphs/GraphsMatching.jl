"""
    minimum_weight_perfect_matching(g, w::Dict{Edge,U<:Real})
    minimum_weight_perfect_matching(g, w::Dict{Edge,U<:Real}, cutoff)

Given a graph `g` and an edgemap `w` containing weights associated to edges,
returns a matching with the mimimum total weight among the ones containing
exactly `nv(g)/2` edges.

Edges in `g` not present in `w` will not be considered for the matching.

This function relies on the BlossomV.jl package, a julia wrapper
around Kolmogorov's BlossomV algorithm.

Eventually a `cutoff` argument can be given, to the reduce computational time
excluding edges with weights higher than the cutoff.

The returned object is of type `MatchingResult`.

If there is no matching, the returned weight will be `typemax(U)` and no vertex will be matched

In case of error try to change the optional argument `tmaxscale` (default is `tmaxscale=2*nv(g)`).
"""
function minimum_weight_perfect_matching end

function minimum_weight_perfect_matching(g::Graph, w::Dict{E,U}, cutoff, kws...) where {U<:Real, E<:Edge}
    wnew = Dict{E, U}()
    for (e, c) in w
        if c <= cutoff
            wnew[e] = c
        end
    end
    return minimum_weight_perfect_matching(g, wnew; kws...)
end

function minimum_weight_perfect_matching(g::Graph, w::Dict{E,U}; tmaxscale=2*nv(g)) where {U<:AbstractFloat, E<:Edge}
    wnew = Dict{E, Int32}()
    cmax = maximum(values(w))
    cmin = minimum(values(w))
    
    
    tmax = typemax(Int32)  / ( (cmax-cmin) * tmaxscale) # tmaxscale = 2*nv(g) is made to ensure that the weights are sufficiently small
                                                    # compared to typemax(Int32)/3, which is used as an infinite value for non edges
    for (e, c) in w
        wnew[e] = round(Int32, (c-cmin) * tmax)
    end
    match = minimum_weight_perfect_matching(g, wnew)
    if match.mate[1] == -1 # there is no matching
        return MatchingResult(typemax(U), match.mate)
    end
    weight = zero(U)
    for i=1:nv(g)
        j = match.mate[i]
        if j > i
            weight += w[E(i,j)]
        end
    end
    return MatchingResult(weight, match.mate)
end

function minimum_weight_perfect_matching(g::Graph, w::Dict{E,U}) where {U<:Integer, E<:Edge}
    @assert nv(g) % 2 == 0
    m = BlossomV.Matching(nv(g))
    for (e, c) in w
        BlossomV.add_edge(m, src(e)-1, dst(e)-1, c)
    end

    # Blossom V needs a feasible matching to work, so we add dummy edges
    tmax = round(Int, typemax(Int32) / 3)
    for i in 1:round(Int, nv(g)/2)
        u, v = 2i-1, 2i
        if Edge(u, v) ∉ keys(w) && Edge(v, u) ∉ keys(w)
            BlossomV.add_edge(m, u-1, v-1, tmax)
        end
    end

    BlossomV.solve(m)

    mate = fill(-1, nv(g))
    totweight = zero(U)
    for i=1:nv(g)
        j = BlossomV.get_match(m, i-1) + 1
        if j == -1 # there is no matching, so we return an infinite weight
            return MatchingResult(typemax(U), fill(-1, nv(g)))
        end
        mate[i] = j <= 0 ? -1 : j
        if i < j
            if !haskey(w, Edge(i, j)) # there is no matching, so we return an infinite weight
                return MatchingResult(typemax(U), fill(-1, nv(g)))
            end
            totweight += w[Edge(i,j)]
        end
    end
    return MatchingResult(totweight, mate)
end
