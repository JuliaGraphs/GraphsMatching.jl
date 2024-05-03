"""
    minimum_weight_perfect_matching(g, w::AbstractMatrix{U})

Given a graph `g` and weights `w`, returns a matching with the mimimum total
weight among the ones containing exactly `nv(g)/2` edges.

This function relies on the BlossomV.jl package, a julia wrapper
around Kolmogorov's BlossomV algorithm.

The returned object is of type `MatchingResult`.

In case of error try to change the optional argument `tmaxscale` (default is `tmaxscale=10`).
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

function minimum_weight_perfect_matching(g::Graph, w::Dict{E,U}; tmaxscale=10.) where {U<:AbstractFloat, E<:Edge}
    wnew = Dict{E, Int32}()
    cmax = maximum(values(w))
    cmin = minimum(values(w))
    tmax = typemax(Int32)  / tmaxscale # /10 is kinda arbitrary,
                                # hopefully high enough to not occur in overflow problems
    for (e, c) in w
        wnew[e] = round(Int32, (c-cmin) / (cmax-cmin) * tmax)
    end
    match = minimum_weight_perfect_matching(g, wnew)
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
    m = BlossomV.Matching(nv(g))
    for (e, c) in w
        BlossomV.add_edge(m, src(e)-1, dst(e)-1, c)
    end
    BlossomV.solve(m)

    mate = fill(-1, nv(g))
    totweight = zero(U)
    for i=1:nv(g)
        j = BlossomV.get_match(m, i-1) + 1
        mate[i] = j <= 0 ? -1 : j
        if i < j
            totweight += w[Edge(i,j)]
        end
    end
    return MatchingResult(totweight, mate)
end

function minimum_weight_perfect_matching(g::Graph, w::AbstractMatrix{U}, cutoff, kws...) where {U <: Real}
	wnew = Dict{Edge, U}()
	for e in edges(g)
		c = w[src(e), dst(e)]
		if c <= cutoff
			wnew[e] = c
		end
	end
	return minimum_weight_perfect_matching(g, wnew, kws...)
end

function minimum_weight_perfect_matching(g::Graph, w::AbstractMatrix{U} = default_weights(g), kws...) where {U <: Real}
	wnew = Dict{Edge, U}()
	for e in edges(g)
		wnew[e] = w[src(e), dst(e)]
	end
	return minimum_weight_perfect_matching(g, wnew, kws...)
end
