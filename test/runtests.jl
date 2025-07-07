using Graphs
using GraphsMatching
using Test
using Cbc
using JuMP
using LinearAlgebra: I

@static if !Sys.iswindows() && Sys.ARCH == :x86_64
    using Pkg
    Pkg.add("BlossomV")
    import BlossomV # to test the extension
end

@testset "GraphsMatching" begin



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
            optimizer = optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
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
            optimizer = optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
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
            optimizer = optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
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
            optimizer = optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0),
            cutoff = 0,
        )
        @test match.weight ≈ 12
        @test match.mate[1] == 6
        @test match.mate[2] == 5
        @test match.mate[3] == -1
        @test match.mate[4] == -1
        @test match.mate[5] == 2
        @test match.mate[6] == 1
    
    end

end
