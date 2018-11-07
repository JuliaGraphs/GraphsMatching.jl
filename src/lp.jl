function maximum_weight_maximal_matching_lp(g::Graph, solver::AbstractMathProgSolver, w::AbstractMatrix{T}, cutoff::R) where {T<:Real, R<:Real}
    return maximum_weight_maximal_matching_lp(g, solver, cutoff_weights(w, cutoff))
end

function maximum_weight_maximal_matching_lp(g::Graph, solver::AbstractMathProgSolver, w::AbstractMatrix{T}) where {T<:Real}
# TODO support for graphs with zero degree nodes
# TODO apply separately on each connected component
    bpmap = bipartite_map(g)
    length(bpmap) != nv(g) && error("Graph is not bipartite")
    v1 = findall(isequal(1), bpmap)
    v2 = findall(isequal(2), bpmap)
    if length(v1) > length(v2)
        v1, v2 = v2, v1
    end

    nedg = 0
    edgemap = Dict{Edge,Int}()

    for j in 1:size(w,2)
        for i in 1:size(w,1)
            if w[i,j] > 0.0
                nedg += 1
                edgemap[Edge(i,j)] = nedg
                edgemap[Edge(j,i)] = nedg
            end
        end
    end

    model = Model(solver=solver)
    @variable(model, x[1:length(w)] >= 0)

    for i in v1
        idx = Vector{Int}()
        for j in neighbors(g, i)
            if haskey(edgemap, Edge(i,j))
                push!(idx, edgemap[Edge(i,j)])
            end
        end
        if length(idx) > 0
            @constraint(model, sum(x[id] for id = idx) == 1)
        end
    end

    for j in v2
        idx = Vector{Int}()
        for i in neighbors(g, j)
            if haskey(edgemap, Edge(i,j))
                push!(idx, edgemap[Edge(i,j)])
            end
        end

        if length(idx) > 0
            @constraint(model, sum(x[id] for id = idx) <= 1)
        end
    end

    @objective(model, Max, sum(w[src(e),dst(e)] * x[edgemap[e]] for e in keys(edgemap)))

    status = solve(model)
    status != :Optimal && error("JuMP solver failed to find optimal solution.")
    sol = getvalue(x)

    all(Bool[s == 1 || s == 0 for s in sol]) || error("Found non-integer solution.")

    cost = getobjectivevalue(model)

    mate = fill(-1, nv(g))
    for e in edges(g)
        if w[src(e),dst(e)] > zero(T)
            inmatch = convert(Bool, sol[edgemap[e]])
            if inmatch
                mate[src(e)] = dst(e)
                mate[dst(e)] = src(e)
            end
        end
    end

    return MatchingResult(cost, mate)
end
