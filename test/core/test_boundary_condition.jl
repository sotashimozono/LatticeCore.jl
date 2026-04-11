using LatticeCore
using Test

@testset "BoundaryCondition basics" begin
    @test PBC <: AbstractBoundaryCondition
    @test OBC <: AbstractBoundaryCondition
    @test PBC() isa AbstractBoundaryCondition
    @test OBC() isa AbstractBoundaryCondition
    @test PBC() != OBC()
end
