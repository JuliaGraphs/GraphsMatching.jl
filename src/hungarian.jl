function maximum_weight_maximal_matching_hungarian(g::Graph,
          w::AbstractMatrix{T}=default_weights(g)) where {T <: Real}
  edge_list = collect(edges(g))
  n = nv(g)

  # Determine the bipartition of the graph. 
  bipartition = bipartite_map(g)
  if length(bipartition) != n # Equivalent to !is_bipartite(g), but reuses the results of the previous function call.
    error("The Hungarian algorithm only works for bipartite graphs; otherwise, prefer the Blossom algorithm (not yet available in LightGraphsMatching")
  end
  n_first = count(bipartition .== 1)
  n_second = count(bipartition .== 2)

  to_bipartition_1 = [count(bipartition[1:i] .== 1) for i in 1:n]
  to_bipartition_2 = [count(bipartition[1:i] .== 2) for i in 1:n]

  # hungarian() minimises the total cost, while this function is supposed to maximise the total weights.
  wDual = maximum(w) .- w

  # Remove weights that are not in the graph (Hungarian.jl considers all weights that are not missing values as real edges).
  # Assume w is symmetric, so that the weight of matching i->j is the same as the one for j->i. 
  weights = Matrix{Union{Missing, T}}(missing, n_first, n_second)

  for i in 1:n
    for j in 1:n
      if Edge(i, j) ∈ edge_list || Edge(j, i) ∈ edge_list 
        if bipartition[i] == 1 # and bipartition[j] == 2
          idx_first = to_bipartition_1[i]
          idx_second = to_bipartition_2[j]
        else # bipartition[i] == 2 and bipartition[j] == 1
          idx_first = to_bipartition_1[j]
          idx_second = to_bipartition_2[i]
        end

        weight_to_add = (Edge(i, j) ∈ edge_list) ? wDual[i, j] : wDual[j, i]

        weights[idx_first, idx_second] = weight_to_add
      end
    end
  end

  # Run the Hungarian algorithm.
  assignment, _ = hungarian(weights)

  # Convert the output format to match LGMatching's. 
  pairs = Tuple{Int, Int}[]
  mate = fill(-1, n) # Initialise to unmatched. 
  for i in eachindex(assignment)
    if assignment[i] != 0 # If matched: 
      original_i = findfirst(to_bipartition_1 .== i)
      original_j = findfirst(to_bipartition_2 .== assignment[i])

      mate[original_i] = original_j
      mate[original_j] = original_i

      push!(pairs, (original_i, original_j))
    end
  end

  # Compute the cost for this matching (as weights had to be changed for Hungarian.jl, the one returned by hungarian() makes no sense). 
  cost = sum(w[p[1], p[2]] for p in pairs)

  # Return the result.
  return MatchingResult(cost, mate)
end
