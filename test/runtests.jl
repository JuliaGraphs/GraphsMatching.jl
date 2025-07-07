using Graphs
using GraphsMatching
using Test
using HiGHS
using JuMP
using LinearAlgebra: I

if !Sys.iswindows() && Sys.ARCH == :x86_64
    using Pkg
    Pkg.add("BlossomV")
    import BlossomV # to test the extension
end

@testset "GraphsMatching" begin

    @testset "maximum_weight_matching_reduction" begin
        N = Int(1e2)
        g = complete_graph(N)
        for i = 1:10
            w = randn(N, N)
            w = transpose(w) * w
            match_sol = maximum_weight_matching(
                g,
                with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
                w,
            )
            match_red = maximum_weight_matching_reduction(g, w)
            sol = Edge[]
            for i = 1:N
                if (match_sol.mate[i] > 0)
                    push!(sol, Edge(i, match_sol.mate[i]))
                end
            end
            @test sol == match_red
        end
    end

    @testset "maximum_weight_matching" begin
        g = complete_graph(3)
        w = [
            1 2 1
            1 1 1
            3 1 1
        ]
        match = maximum_weight_matching(
            g,
            with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            w,
        )
        @test match.mate[1] == 3
        @test match.weight ≈ 3

        g = complete_graph(3)
        w = zeros(3, 3)
        w[1, 2] = 1
        w[3, 2] = 1
        w[1, 3] = 1
        match = maximum_weight_matching(
            g,
            with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            w,
        )
        @test match.weight ≈ 1


        g = Graph(4)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 4)

        w = zeros(4, 4)
        w[1, 3] = 1
        w[1, 4] = 3
        w[2, 4] = 1

        match = maximum_weight_matching(
            g,
            with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            w,
        )
        @test match.weight ≈ 3
        @test match.mate[1] == 4
        @test match.mate[2] == -1
        @test match.mate[3] == -1
        @test match.mate[4] == 1

        g = Graph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        add_edge!(g, 3, 4)
        match = maximum_weight_matching(
            g,
            with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
        )
        @test match.weight ≈ 2
        @test match.mate[1] == 2
        @test match.mate[2] == 1
        @test match.mate[3] == 4
        @test match.mate[4] == 3

        w = zeros(4, 4)
        w[1, 2] = 1
        w[2, 3] = 1
        w[1, 3] = 1
        w[3, 4] = 1

        match = maximum_weight_matching(
            g,
            with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            w,
        )
        @test match.weight ≈ 2
        @test match.mate[1] == 2
        @test match.mate[2] == 1
        @test match.mate[3] == 4
        @test match.mate[4] == 3

        w = zeros(4, 4)
        w[1, 2] = 1
        w[2, 3] = 1
        w[1, 3] = 5
        w[3, 4] = 1

        match = maximum_weight_matching(
            g,
            with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            w,
        )
        @test match.weight ≈ 5
        @test match.mate[1] == 3
        @test match.mate[2] == -1
        @test match.mate[3] == 1
        @test match.mate[4] == -1
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
            w,
            algorithm = LPAlgorithm(),
            optimizer = with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
        )
        @test match.weight ≈ 21
        @test match.mate[1] == 3
        @test match.mate[3] == 1
        @test match.mate[2] == 4
        @test match.mate[4] == 2

        g = complete_bipartite_graph(2, 4)
        w = zeros(6, 6)
        w[1, 3] = 10
        w[1, 4] = 0.5
        w[2, 3] = 11
        w[2, 4] = 1
        match = maximum_weight_maximal_matching(
            g,
            w,
            algorithm = LPAlgorithm(),
            optimizer = with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
        )
        @test match.weight ≈ 11.5
        @test match.mate[1] == 4
        @test match.mate[4] == 1
        @test match.mate[2] == 3
        @test match.mate[3] == 2

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
            w,
            algorithm = LPAlgorithm(),
            optimizer = with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            cutoff = 0,
        )
        @test match.weight ≈ 11.5
        @test match.mate[1] == 4
        @test match.mate[4] == 1
        @test match.mate[2] == 3
        @test match.mate[3] == 2

        g = complete_bipartite_graph(4, 2)
        w = zeros(6, 6)
        w[3, 5] = 10
        w[3, 6] = 0.5
        w[2, 5] = 11
        w[1, 6] = 1
        w[1, 5] = -1

        match = maximum_weight_maximal_matching(
            g,
            w,
            algorithm = LPAlgorithm(),
            optimizer = with_optimizer(HiGHS.Optimizer, "logLevel" => 0),
            cutoff = 0,
        )
        @test match.weight ≈ 12
        @test match.mate[1] == 6
        @test match.mate[2] == 5
        @test match.mate[3] == -1
        @test match.mate[4] == -1
        @test match.mate[5] == 2
        @test match.mate[6] == 1


        g = complete_bipartite_graph(2, 2)
        w = zeros(4, 4)
        w[1, 3] = 10.0
        w[1, 4] = 1.0
        w[2, 3] = 2.0
        w[2, 4] = 11.0
        match = maximum_weight_maximal_matching(g, w, algorithm = HungarianAlgorithm())
        @test match.weight ≈ 21
        @test match.mate[1] == 3
        @test match.mate[3] == 1
        @test match.mate[2] == 4
        @test match.mate[4] == 2

        g = complete_graph(3)
        w = zeros(3, 3)
        @test !is_bipartite(g)
        @test_throws ErrorException maximum_weight_maximal_matching(
            g,
            w,
            algorithm = HungarianAlgorithm(),
        )

        g = complete_bipartite_graph(2, 4)
        w = zeros(6, 6)
        w[1, 3] = 10
        w[1, 4] = 0.5
        w[2, 3] = 11
        w[2, 4] = 1
        match = maximum_weight_maximal_matching(g, w, algorithm = HungarianAlgorithm())
        @test match.weight ≈ 11.5

        g = Graph(4)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        add_edge!(g, 2, 4)
        w = zeros(4, 4)
        w[1, 3] = 1
        w[1, 4] = 3
        w[2, 4] = 1
        match = maximum_weight_maximal_matching(g, w, algorithm = HungarianAlgorithm())
        @test match.weight ≈ 2

    end



    @testset "minimum_weight_perfect_matching" begin
        algos = Any[LEMONMWPMAlgorithm()]
        if !Sys.iswindows() && Sys.ARCH == :x86_64
            push!(algos, BlossomVAlgorithm())
        end
        for algorithm in algos

        w = Dict(Edge(1, 2) => 500)
        g = Graph(2)
        add_edge!(g, 1, 2)
        match = minimum_weight_perfect_matching(g, w, algorithm)
        @test match.mate[1] == 2


        w = Dict(
            Edge(1, 2) => 500,
            Edge(1, 3) => 600,
            Edge(2, 3) => 700,
            Edge(3, 4) => 100,
            Edge(2, 4) => 1000,
        )

        g = complete_graph(4)
        match = minimum_weight_perfect_matching(g, w, algorithm)
        @test match.mate[1] == 2
        @test match.mate[2] == 1
        @test match.mate[3] == 4
        @test match.mate[4] == 3
        @test match.weight ≈ 600

        w = Dict(
            Edge(1, 2) => 500,
            Edge(1, 3) => 400,
            Edge(2, 3) => 300,
            Edge(3, 4) => 1000,
            Edge(2, 4) => 1000,
        )
        g = complete_graph(4)
        match = minimum_weight_perfect_matching(g, w, algorithm)
        @test match.mate[1] == 3
        @test match.mate[2] == 4
        @test match.mate[3] == 1
        @test match.mate[4] == 2
        @test match.weight ≈ 1400

        g = complete_bipartite_graph(2, 2)
        w = Dict{Edge,Float64}()
        w[Edge(1, 3)] = -10
        w[Edge(1, 4)] = -0.5
        w[Edge(2, 3)] = -11
        w[Edge(2, 4)] = -1

        match = minimum_weight_perfect_matching(g, w, algorithm)
        @test match.mate[1] == 4
        @test match.mate[4] == 1
        @test match.mate[2] == 3
        @test match.mate[3] == 2
        @test match.weight ≈ -11.5


        g = complete_graph(4)
        w = Dict{Edge,Float64}()
        w[Edge(1, 3)] = 10
        w[Edge(1, 4)] = 0.5
        w[Edge(2, 3)] = 11
        w[Edge(2, 4)] = 2
        w[Edge(1, 2)] = 100

        match = minimum_weight_perfect_matching(g, w, 50, algorithm)
        @test match.mate[1] == 4
        @test match.mate[4] == 1
        @test match.mate[2] == 3
        @test match.mate[3] == 2
        @test match.weight ≈ 11.5

        # Issue #5
        g = Graph(Graph([Edge(1, 2), Edge(2, 3), Edge(4, 1)]))
        w = Dict(Edge(1, 2) => 2.0, Edge(2, 3) => 2.0, Edge(1, 4) => 2.0)
        match = minimum_weight_perfect_matching(g, w)
        @test match.mate == [4, 3, 2, 1]
        @test match.weight == 4.0

        g = Graph([Edge(1, 2)])
        wFloat = Dict(Edge(1, 2) => 2.0)
        wInt = Dict(Edge(1, 2) => 2)
        matchFloat = minimum_weight_perfect_matching(g, wFloat, algorithm)
        matchInt = minimum_weight_perfect_matching(g, wInt, algorithm)
        @test matchFloat.mate == matchInt.mate
        @test matchFloat.weight == matchInt.weight

        end
    end
end
