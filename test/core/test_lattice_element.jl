using LatticeCore
using StaticArrays
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

@testset "Element-center generic API on SimpleSquareLattice" begin
    lat = SimpleSquareLattice(3, 3, OpenAxis())
    nbond = count(_ -> true, bonds(lat))

    @testset "num_elements" begin
        @test num_elements(lat, VertexCenter()) == num_sites(lat)
        @test num_elements(lat, BondCenter()) == nbond
        @test_throws MethodError num_elements(lat, PlaquetteCenter())
    end

    @testset "elements iterator" begin
        verts = collect(elements(lat, VertexCenter()))
        @test verts == 1:num_sites(lat)

        bs = collect(elements(lat, BondCenter()))
        @test length(bs) == nbond
        @test all(b isa Bond{2,Float64} for b in bs)
    end

    @testset "element_position dispatches by centring" begin
        @test element_position(lat, VertexCenter(), 1) == position(lat, 1)

        bs = collect(bonds(lat))
        @test element_position(lat, BondCenter(), 1) == bond_center(lat, bs[1])
    end

    @testset "element_positions iterator" begin
        vps = collect(element_positions(lat, VertexCenter()))
        @test vps == [position(lat, i) for i in 1:num_sites(lat)]

        bps = collect(element_positions(lat, BondCenter()))
        @test length(bps) == nbond
        @test all(p isa SVector{2,Float64} for p in bps)
    end
end
