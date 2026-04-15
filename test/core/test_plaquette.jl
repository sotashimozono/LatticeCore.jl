using LatticeCore
using StaticArrays
using Test

# Minimal concrete lattice used to exercise the plaquette /
# element_neighbors / incident defaults without depending on
# downstream packages. Models a 2×2 square sample with OBC so we can
# count everything by hand.
struct TinySquareLattice <: AbstractLattice{2,Float64} end

LatticeCore.num_sites(::TinySquareLattice) = 4
function LatticeCore.position(::TinySquareLattice, i::Int)
    (
        if (i == 1)
            SVector(0.0, 0.0)
        elseif (i == 2)
            SVector(1.0, 0.0)
        elseif (i == 3)
            SVector(0.0, 1.0)
        else
            SVector(1.0, 1.0)
        end
    )
end
function LatticeCore.neighbors(::TinySquareLattice, i::Int)
    if i == 1
        [2, 3]
    elseif i == 2
        [1, 4]
    elseif i == 3
        [1, 4]
    else
        [2, 3]
    end
end
LatticeCore.boundary(::TinySquareLattice) = nothing
LatticeCore.size_trait(::TinySquareLattice) = FiniteSize((2, 2))

# One square plaquette covering all four sites, cyclic order 1→2→4→3.
function LatticeCore.plaquettes(::TinySquareLattice)
    (Plaquette{2,Float64}([1, 2, 4, 3], SVector(0.5, 0.5), :square),)
end

@testset "PlaquetteRule and Plaquette types" begin
    rule = PlaquetteRule([(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)], :square)
    @test rule.corners[1] == (1, 0, 0)
    @test rule.type === :square
    @test length(rule.corners) == 4

    p = Plaquette{2,Float64}([1, 2, 4, 3], SVector(0.5, 0.5), :square)
    @test p.vertices == [1, 2, 4, 3]
    @test plaquette_center(p) == SVector(0.5, 0.5)
    @test p.type === :square
end

@testset "Element-center defaults for PlaquetteCenter" begin
    lat = TinySquareLattice()

    @test num_elements(lat, PlaquetteCenter()) == 1
    @test length(collect(elements(lat, PlaquetteCenter()))) == 1
    @test element_position(lat, PlaquetteCenter(), 1) == SVector(0.5, 0.5)
end

@testset "neighbor_plaquettes default" begin
    lat = TinySquareLattice()
    # Every site is on the boundary of the single square plaquette.
    for i in 1:4
        @test length(collect(neighbor_plaquettes(lat, i))) == 1
    end
end

@testset "element_neighbors (line graph / dual graph)" begin
    lat = TinySquareLattice()

    @testset "VertexCenter → existing neighbors" begin
        @test element_neighbors(lat, VertexCenter(), 1) == neighbors(lat, 1)
    end

    @testset "BondCenter → line graph: bonds sharing a vertex" begin
        bs = collect(bonds(lat))
        # TinySquare has 4 bonds: 1-2, 1-3, 2-4, 3-4.
        @test length(bs) == 4
        # Every bond shares a vertex with exactly two other bonds in
        # this tiny square.
        for i in 1:length(bs)
            @test length(element_neighbors(lat, BondCenter(), i)) == 2
        end
    end

    @testset "PlaquetteCenter → dual graph: single plaquette has no peers" begin
        @test isempty(element_neighbors(lat, PlaquetteCenter(), 1))
    end
end

@testset "incident: cross-centring" begin
    lat = TinySquareLattice()
    bs = collect(bonds(lat))

    @testset "VertexCenter → BondCenter: site i is on exactly 2 bonds" begin
        for i in 1:4
            inc = incident(lat, VertexCenter(), BondCenter(), i)
            @test length(inc) == 2
            # Every returned bond index must actually touch site i.
            @test all(b -> b.i == i || b.j == i, bs[inc])
        end
    end

    @testset "BondCenter → VertexCenter: endpoints" begin
        for k in 1:length(bs)
            @test Set(incident(lat, BondCenter(), VertexCenter(), k)) ==
                Set([bs[k].i, bs[k].j])
        end
    end

    @testset "VertexCenter → PlaquetteCenter: every site is on the plaquette" begin
        for i in 1:4
            @test incident(lat, VertexCenter(), PlaquetteCenter(), i) == [1]
        end
    end

    @testset "PlaquetteCenter → VertexCenter: boundary walk" begin
        @test incident(lat, PlaquetteCenter(), VertexCenter(), 1) == [1, 2, 4, 3]
    end

    @testset "BondCenter → PlaquetteCenter: every bond bounds the plaquette" begin
        for k in 1:length(bs)
            @test incident(lat, BondCenter(), PlaquetteCenter(), k) == [1]
        end
    end

    @testset "PlaquetteCenter → BondCenter: the 4 boundary bonds" begin
        @test length(incident(lat, PlaquetteCenter(), BondCenter(), 1)) == 4
    end

    @testset "Same-centring pairs fall through to element_neighbors" begin
        @test incident(lat, VertexCenter(), VertexCenter(), 1) == neighbors(lat, 1)
    end
end
