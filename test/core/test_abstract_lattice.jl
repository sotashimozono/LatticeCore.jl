using LatticeCore
using StaticArrays
using Test

# A minimal concrete subtype used to exercise the abstract interface.
# Deliberately local to this test file: the "real" reference
# implementations (LineLattice, SimpleSquareLattice) will land in a
# follow-up PR.

struct TrivialBoundary end

struct TestLattice1D <: AbstractLattice{1,Float64}
    N::Int
end

LatticeCore.num_sites(l::TestLattice1D) = l.N
LatticeCore.position(l::TestLattice1D, i::Int) = SVector{1,Float64}(Float64(i))
function LatticeCore.neighbors(l::TestLattice1D, i::Int)
    return filter(j -> 1 <= j <= l.N, [i - 1, i + 1])
end
LatticeCore.boundary(::TestLattice1D) = TrivialBoundary()
LatticeCore.size_trait(l::TestLattice1D) = FiniteSize((l.N,))

@testset "AbstractLattice via TestLattice1D" begin
    lat = TestLattice1D(5)

    @testset "type parameters and basic accessors" begin
        @test lat isa AbstractLattice{1,Float64}
        @test dimension(lat) == 1
        @test eltype(lat) === Float64
        @test eltype(typeof(lat)) === Float64
        @test num_sites(lat) == 5
    end

    @testset "position" begin
        @test position(lat, 1) == SVector(1.0)
        @test position(lat, 3) == SVector(3.0)
        @test position(lat, 5) == SVector(5.0)
    end

    @testset "neighbors" begin
        @test neighbors(lat, 1) == [2]
        @test neighbors(lat, 3) == [2, 4]
        @test neighbors(lat, 5) == [4]
    end

    @testset "boundary" begin
        @test boundary(lat) isa TrivialBoundary
    end

    @testset "size trait defaults and is_finite" begin
        st = size_trait(lat)
        @test st isa FiniteSize{1}
        @test st.dims == (5,)
        @test is_finite(lat) == true
    end

    @testset "trait defaults on AbstractLattice" begin
        @test topology(lat) isa TopologyTrait{:unknown}
        @test periodicity(lat) isa Aperiodic
        @test is_bipartite(lat) == false
        @test reciprocal_support(lat) isa NoReciprocal
    end

    @testset "positions iterator (default)" begin
        ps = collect(positions(lat))
        @test length(ps) == 5
        @test ps[1] == SVector(1.0)
        @test ps[5] == SVector(5.0)
    end

    @testset "bonds iterator (default)" begin
        all_bonds = collect(bonds(lat))
        # Expected: (1,2), (2,3), (3,4), (4,5) — 4 bonds for a 1D chain
        @test length(all_bonds) == 4
        @test all(b isa Bond{1,Float64} for b in all_bonds)
        @test all(b.type === :nearest for b in all_bonds)
        @test all_bonds[1].i == 1 && all_bonds[1].j == 2
        @test all_bonds[end].i == 4 && all_bonds[end].j == 5
    end

    @testset "neighbor_bonds iterator (default)" begin
        nb_mid = collect(neighbor_bonds(lat, 3))
        @test length(nb_mid) == 2
        @test Set((b.i, b.j) for b in nb_mid) == Set([(3, 2), (3, 4)])

        nb_left = collect(neighbor_bonds(lat, 1))
        @test length(nb_left) == 1
        @test nb_left[1].j == 2
    end

    @testset "bond_center" begin
        b = Bond{1,Float64}(1, 2, SVector(1.0), :nearest)
        @test bond_center(lat, b) == SVector(1.5)

        b2 = Bond{1,Float64}(3, 4, SVector(1.0), :nearest)
        @test bond_center(lat, b2) == SVector(3.5)
    end
end
