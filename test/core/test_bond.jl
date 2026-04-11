using LatticeCore
using StaticArrays
using Test

@testset "Bond struct" begin
    b = Bond{2, Float64}(1, 2, SVector(1.0, 0.0), :nearest)
    @test b.i == 1
    @test b.j == 2
    @test b.vector == SVector(1.0, 0.0)
    @test b.type === :nearest

    # Different bond type tag
    b2 = Bond{2, Float64}(3, 7, SVector(-0.5, 0.5), :dimer_strong)
    @test b2.type === :dimer_strong
    @test b2.vector == SVector(-0.5, 0.5)

    # Bond is immutable (typeof check)
    @test isimmutable(b)
end
