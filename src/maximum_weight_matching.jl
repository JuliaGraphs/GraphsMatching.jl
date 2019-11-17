"""
maximum_weight_matching(g::Graph, w::Dict{Edge,Real} -> Dict{Edge,Int64}

Given a graph `g` and an edgemap `w` containing weights associated to edges,
returns a matching with the maximum total weight.
`w` is a dictionary that maps edges i => j to weights.
If no weight parameter is given, all edges will be considered to have weight 1
(results in max cardinality matching)

The efficiency of the algorithm depends on the input graph:
  - If the graph is bipartite, then the LP relaxation is integral.
  - If the graph is not bipartite, then it requires a MIP solver and
  the computation time may grow exponentially.

The package JuMP.jl and one of its supported solvers is required.

Returns MatchingResult containing:
  - a solve status (indicating whether the problem was solved to optimality)
  - the optimal cost
  - a list of each vertex's match (or -1 for unmatched vertices)
"""
function maximum_weight_matching end

function maximum_weight_matching(g::Graph,
          solver::JuMP.OptimizerFactory,
          w::AbstractMatrix{U} = default_weights(g)) where {U <:Real}

    model = Model(with_optimizer(solver))
    n = nv(g)
    edge_list = collect(edges(g))

    # put the edge weights in w in the right order to be compatible with edge_list
    for j in 1:n
      for i in 1:n
        if i > j  && w[i,j] > zero(U) && w[j,i] < w[i,j]
          w[j,i] = w[i,j]
        end
        if Edge(i,j) ∉ edge_list
          w[i,j] = zero(U)
        end
      end
    end
    
    if is_bipartite(g)
      @variable(model, x[edge_list] >= 0) # no need to enforce integrality
    else
      @variable(model, x[edge_list] >= 0, Int) # requires MIP solver
    end
    @objective(model, Max, sum(x[e]*w[src(e),dst(e)] for e in edge_list))

    @constraint(model, c1[i=1:n], sum(x[Edge(minmax(i,j))] for j in neighbors(g,i)) <= 1)
    optimize!(model)
    status = JuMP.termination_status(model)
    status != MOI.OPTIMAL && error("JuMP solver failed to find optimal solution.")
    solution = value.(x)
    cost = objective_value(model)
    return MatchingResult(cost, dict_to_arr(n, solution, edge_list))
end

""" Returns an array of mates from a dictionary that maps edges to {0,1} """
function dict_to_arr(n::Int64, solution::JuMP.Containers.DenseAxisArray{U,1,Tuple{Array{E,1}}}, edge_list::AbstractVector{E}) where {U<: Real, E<: Edge}
  mate = fill(-1,n)
  for e in edge_list
    if solution[e] >= 1 - 1e-5 # Some tolerance to numerical approximations by the solver.
        mate[src(e)] = dst(e)
        mate[dst(e)] = src(e)
    end
  end
  return mate
end


function default_weights(g::G) where {G<:AbstractGraph}
  m = spzeros(nv(g),nv(g))
  for e in edges(g)
    m[src(e),dst(e)] = 1
  end
  return m
end
