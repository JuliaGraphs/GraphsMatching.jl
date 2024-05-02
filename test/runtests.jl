using Graphs
using GraphsMatching
using Test
using Cbc
using JuMP
using LinearAlgebra: I

@testset "GraphsMatching" begin
    @testset "accessors" begin
        match = MatchingResult(45, [4, -1, 7, 1, 6, 5, 3, -1])

        @test weight(match) == 45

        @test is_matched_vertex(match, 3) == true
        @test is_matched_vertex(match, 2) == false

        @test matching_vertex(match, 5) == 6
        @test matching_vertex(match, 8) == -1

        @test matching_vertices(match) == [4, -1, 7, 1, 6, 5, 3, -1]

        m_edges = matched_edges(match)
        @test length(m_edges) == 3
        @test Edge(1, 4) ∈ m_edges
        @test Edge(3, 7) ∈ m_edges
        @test Edge(5, 6) ∈ m_edges
    end

    @testset "maximum_weight_matching" begin
        g = complete_graph(3)
        w = [
            1 2 1
            1 1 1
            3 1 1
        ]
        match = maximum_weight_matching(
            g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0), w
        )
        @test matching_vertex(match, 1) == 3
        @test weight(match) ≈ 3

        g = complete_graph(3)
        w = zeros(3, 3)
        w[1, 2] = 1
        w[3, 2] = 1
        w[1, 3] = 1
        match = maximum_weight_matching(
            g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0), w
        )
        @test weight(match) ≈ 1

        g = Graph(4)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 4)

        w = zeros(4, 4)
        w[1, 3] = 1
        w[1, 4] = 3
        w[2, 4] = 1

        match = maximum_weight_matching(
            g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0), w
        )
        @test weight(match) ≈ 3
        @test matching_vertex(match, 1) == 4
        @test matching_vertex(match, 2) == -1
        @test matching_vertex(match, 3) == -1
        @test matching_vertex(match, 4) == 1

        g = Graph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        add_edge!(g, 3, 4)
        match = maximum_weight_matching(
            g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0)
        )
        @test weight(match) ≈ 2
        @test matching_vertex(match, 1) == 2
        @test matching_vertex(match, 2) == 1
        @test matching_vertex(match, 3) == 4
        @test matching_vertex(match, 4) == 3

        w = zeros(4, 4)
        w[1, 2] = 1
        w[2, 3] = 1
        w[1, 3] = 1
        w[3, 4] = 1

        match = maximum_weight_matching(
            g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0), w
        )
        @test weight(match) ≈ 2
        @test matching_vertex(match, 1) == 2
        @test matching_vertex(match, 2) == 1
        @test matching_vertex(match, 3) == 4
        @test matching_vertex(match, 4) == 3

        w = zeros(4, 4)
        w[1, 2] = 1
        w[2, 3] = 1
        w[1, 3] = 5
        w[3, 4] = 1

        match = maximum_weight_matching(
            g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0), w
        )
        @test weight(match) ≈ 5
        @test matching_vertex(match, 1) == 3
        @test matching_vertex(match, 2) == -1
        @test matching_vertex(match, 3) == 1
        @test matching_vertex(match, 4) == -1
    end

    @testset "maximum_weight_maximal_matching" begin
        g = complete_bipartite_graph(2, 2)
        w = zeros(4, 4)
        w[1, 3] = 10.0
        w[1, 4] = 1.0
        w[2, 3] = 2.0
        w[2, 4] = 11.0
        match = maximum_weight_maximal_matching(
            g,
            w;
            algorithm=LPAlgorithm(),
            optimizer=optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
        )
        @test weight(match) ≈ 21
        @test matching_vertex(match, 1) == 3
        @test matching_vertex(match, 3) == 1
        @test matching_vertex(match, 2) == 4
        @test matching_vertex(match, 4) == 2

        g = complete_bipartite_graph(2, 4)
        w = zeros(6, 6)
        w[1, 3] = 10
        w[1, 4] = 0.5
        w[2, 3] = 11
        w[2, 4] = 1
        match = maximum_weight_maximal_matching(
            g,
            w;
            algorithm=LPAlgorithm(),
            optimizer=optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
        )
        @test weight(match) ≈ 11.5
        @test matching_vertex(match, 1) == 4
        @test matching_vertex(match, 4) == 1
        @test matching_vertex(match, 2) == 3
        @test matching_vertex(match, 3) == 2

        g = complete_bipartite_graph(2, 6)
        w = zeros(8, 8)
        w[1, 3] = 10
        w[1, 4] = 0.5
        w[2, 3] = 11
        w[2, 4] = 1
        w[2, 5] = -1
        w[2, 6] = -1
        match = maximum_weight_maximal_matching(
            g,
            w;
            algorithm=LPAlgorithm(),
            optimizer=optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
            cutoff=0,
        )
        @test weight(match) ≈ 11.5
        @test matching_vertex(match, 1) == 4
        @test matching_vertex(match, 4) == 1
        @test matching_vertex(match, 2) == 3
        @test matching_vertex(match, 3) == 2

        g = complete_bipartite_graph(4, 2)
        w = zeros(6, 6)
        w[3, 5] = 10
        w[3, 6] = 0.5
        w[2, 5] = 11
        w[1, 6] = 1
        w[1, 5] = -1

        match = maximum_weight_maximal_matching(
            g,
            w;
            algorithm=LPAlgorithm(),
            optimizer=optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
            cutoff=0,
        )
        @test weight(match) ≈ 12
        @test matching_vertex(match, 1) == 6
        @test matching_vertex(match, 2) == 5
        @test matching_vertex(match, 3) == -1
        @test matching_vertex(match, 4) == -1
        @test matching_vertex(match, 5) == 2
        @test matching_vertex(match, 6) == 1

        g = complete_bipartite_graph(2, 2)
        w = zeros(4, 4)
        w[1, 3] = 10.0
        w[1, 4] = 1.0
        w[2, 3] = 2.0
        w[2, 4] = 11.0
        match = maximum_weight_maximal_matching(g, w; algorithm=HungarianAlgorithm())
        @test weight(match) ≈ 21
        @test matching_vertex(match, 1) == 3
        @test matching_vertex(match, 3) == 1
        @test matching_vertex(match, 2) == 4
        @test matching_vertex(match, 4) == 2

        g = complete_graph(3)
        w = zeros(3, 3)
        @test !is_bipartite(g)
        @test_throws ErrorException maximum_weight_maximal_matching(
            g, w, algorithm=HungarianAlgorithm()
        )

        g = complete_bipartite_graph(2, 4)
        w = zeros(6, 6)
        w[1, 3] = 10
        w[1, 4] = 0.5
        w[2, 3] = 11
        w[2, 4] = 1
        match = maximum_weight_maximal_matching(g, w; algorithm=HungarianAlgorithm())
        @test weight(match) ≈ 11.5

        g = Graph(4)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 4)
        w = zeros(4, 4)
        w[1, 3] = 1
        w[1, 4] = 3
        w[2, 4] = 1
        match = maximum_weight_maximal_matching(g, w; algorithm=HungarianAlgorithm())
        @test weight(match) ≈ 2
    end

    @testset "minimum_weight_perfect_matching" begin
        w = Dict(Edge(1, 2) => 500)
        g = Graph(2)
        add_edge!(g, 1, 2)
        match = minimum_weight_perfect_matching(g, w)
        @test matching_vertex(match, 1) == 2

        w = Dict(
            Edge(1, 2) => 500,
            Edge(1, 3) => 600,
            Edge(2, 3) => 700,
            Edge(3, 4) => 100,
            Edge(2, 4) => 1000,
        )

        g = complete_graph(4)
        match = minimum_weight_perfect_matching(g, w)
        @test matching_vertex(match, 1) == 2
        @test matching_vertex(match, 2) == 1
        @test matching_vertex(match, 3) == 4
        @test matching_vertex(match, 4) == 3
        @test weight(match) ≈ 600

        w = Dict(
            Edge(1, 2) => 500,
            Edge(1, 3) => 400,
            Edge(2, 3) => 300,
            Edge(3, 4) => 1000,
            Edge(2, 4) => 1000,
        )
        g = complete_graph(4)
        match = minimum_weight_perfect_matching(g, w)
        @test matching_vertex(match, 1) == 3
        @test matching_vertex(match, 2) == 4
        @test matching_vertex(match, 3) == 1
        @test matching_vertex(match, 4) == 2
        @test weight(match) ≈ 1400

        g = complete_bipartite_graph(2, 2)
        w = Dict{Edge,Float64}()
        w[Edge(1, 3)] = -10
        w[Edge(1, 4)] = -0.5
        w[Edge(2, 3)] = -11
        w[Edge(2, 4)] = -1

        match = minimum_weight_perfect_matching(g, w)
        @test matching_vertex(match, 1) == 4
        @test matching_vertex(match, 4) == 1
        @test matching_vertex(match, 2) == 3
        @test matching_vertex(match, 3) == 2
        @test weight(match) ≈ -11.5

        g = complete_graph(4)
        w = Dict{Edge,Float64}()
        w[Edge(1, 3)] = 10
        w[Edge(1, 4)] = 0.5
        w[Edge(2, 3)] = 11
        w[Edge(2, 4)] = 2
        w[Edge(1, 2)] = 100

        match = minimum_weight_perfect_matching(g, w, 50)
        @test matching_vertex(match, 1) == 4
        @test matching_vertex(match, 4) == 1
        @test matching_vertex(match, 2) == 3
        @test matching_vertex(match, 3) == 2
        @test weight(match) ≈ 11.5
    end
end
