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
    @test 1==1
end
