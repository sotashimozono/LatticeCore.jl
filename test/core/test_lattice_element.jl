using LatticeCore
using Test

@testset "AbstractLatticeElement" begin
    @test VertexCenter <: AbstractLatticeElement
    @test BondCenter <: AbstractLatticeElement
    @test PlaquetteCenter <: AbstractLatticeElement
    @test CellCenter <: AbstractLatticeElement

    # Singleton values
    @test VertexCenter() isa AbstractLatticeElement
    @test BondCenter() isa AbstractLatticeElement
    @test PlaquetteCenter() isa AbstractLatticeElement
    @test CellCenter() isa AbstractLatticeElement
end
