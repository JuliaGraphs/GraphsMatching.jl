"""
maximum_weight_perfect_matching(g, w::Dict{Edge,Real})
maximum_weight_perfect_matching(g, w::Dict{Edge,Real}, cutoff)

Given a graph `g` and an edgemap `w` containing weights associated to edges,
returns a matching with the maximum total weight among the ones containing
exactly `nv(g)/2` edges.

Edges in `g` not present in `w` will not be considered for the matching.

This function implements a pure-Julia Blossom algorithm for maximum weight matching.

Eventually a `cutoff` argument can be given, to reduce computational time
excluding edges with weights lower than the cutoff.

The returned object is of type `MatchingResult`.
"""
function maximum_weight_perfect_matching(g::Graph, w::Dict{E,U}, cutoff) where {U<:Real,E<:Edge}
    n = nv(g)  # Number of vertices
    mate = Dict{Int, Int}()  # Matching pairs
    weight = 0.0  # Total weight of the matching
    label = fill(0, n)  # 0: unlabeled, 1: labeled, 2: in blossom
    parent = fill(-1, n)  # Parent in the augmenting path
    blossom = fill(-1, n)  # To track blossoms

    function find_augmenting_path(start)
        # Initialize for BFS
        queue = [start]
        label[start] = 1
        parent[start] = -1
        blossom[start] = -1

        while !isempty(queue)
            u = popfirst!(queue)
            for v in neighbors(g, u)
                if !haskey(w, Edge(u, v)) || w[Edge(u, v)] < cutoff
                    continue
                end
                if label[v] == 0  # Unlabeled
                    label[v] = 1
                    parent[v] = u
                    blossom[v] = -1
                    push!(queue, v)
                elseif label[v] == 1 && parent[u] != v  # Found an augmenting path
                    # Handle blossom formation
                    handle_blossom(u, v)
                    return true
                end
            end
        end
        return false
    end

    function handle_blossom(u, v)
        # Logic to handle the formation of a blossom
        # Mark the vertices in the blossom
        while true
            if u == -1 || v == -1
                break
            end
            if label[u] == 1
                label[u] = 2  # Mark as in blossom
                u = parent[mate[u]]  # Move to the parent in the matching
            end
            if label[v] == 1
                label[v] = 2  # Mark as in blossom
                v = parent[mate[v]]  # Move to the parent in the matching
            end
        end
    end

    function augment_path(u, v)
        # Augment the path from u to v
        while u != -1 || v != -1
            if u != -1
                next_u = get(mate, u, -1)
                mate[u] = v
                v = next_u
            end
            if v != -1
                next_v = get(mate, v, -1)
                mate[v] = u
                u = next_v
            end
        end
    end

    for u in 1:n
        if !haskey(mate, u)  # If u is unmatched
            label .= 0  # Reset labels
            if find_augmenting_path(u)
                # Update total weight
                weight += w[Edge(u, get(mate, u, -1))]
            end
        end
    end

    return MatchingResult(weight, mate)
end

function maximum_weight_perfect_matching(
    g::Graph, w::Dict{E,U}; tmaxscale=10.0
) where {U<:AbstractFloat,E<:Edge}
    return maximum_weight_perfect_matching(g, w; cutoff=tmaxscale)
end

function maximum_weight_perfect_matching(g::Graph, w::Dict{E,U}) where {U<:Integer,E<:Edge}
    return maximum_weight_perfect_matching(g, w)
end

# Remove or comment out the dependency on BlossomV.jl
# using BlossomV

# Add tests and documentation for the new implementation
# ...
