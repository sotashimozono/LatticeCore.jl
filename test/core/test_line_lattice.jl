using LatticeCore
using StaticArrays
using Test

@testset "LineLattice" begin
    @testset "construction (default Float64)" begin
        lat = LineLattice(5)
        @test lat isa LineLattice{Float64,PBC}
        @test lat isa AbstractLattice{1,Float64}
        @test dimension(lat) == 1
        @test eltype(lat) === Float64
        @test num_sites(lat) == 5
        @test boundary(lat) isa PBC
    end

    @testset "PBC neighbors" begin
        lat = LineLattice(5, PBC())
        # interior
        @test Set(neighbors(lat, 3)) == Set([2, 4])
        # left edge wraps to right edge
        @test Set(neighbors(lat, 1)) == Set([5, 2])
        # right edge wraps to left edge
        @test Set(neighbors(lat, 5)) == Set([4, 1])
    end

    @testset "OBC neighbors" begin
        lat = LineLattice(5, OBC())
        @test neighbors(lat, 1) == [2]
        @test Set(neighbors(lat, 3)) == Set([2, 4])
        @test neighbors(lat, 5) == [4]
    end

    @testset "positions" begin
        lat = LineLattice(4, OBC())
        ps = collect(positions(lat))
        @test length(ps) == 4
        @test ps[1] == SVector(1.0)
        @test ps[4] == SVector(4.0)
    end

    @testset "bond counts" begin
        # OBC chain of length N has N-1 bonds
        lat_obc = LineLattice(5, OBC())
        @test length(collect(bonds(lat_obc))) == 4
        # PBC cycle of length N has N bonds
        lat_pbc = LineLattice(5, PBC())
        @test length(collect(bonds(lat_pbc))) == 5
    end

    @testset "neighbor_bonds" begin
        lat = LineLattice(5, PBC())
        nb_mid = collect(neighbor_bonds(lat, 3))
        @test length(nb_mid) == 2
        @test all(b.i == 3 for b in nb_mid)
        @test Set(b.j for b in nb_mid) == Set([2, 4])
    end

    @testset "size trait / is_finite" begin
        lat = LineLattice(7, OBC())
        @test size_trait(lat) isa FiniteSize{1}
        @test size_trait(lat).dims == (7,)
        @test is_finite(lat) == true
    end

    @testset "trait overrides" begin
        pbc_even = LineLattice(6, PBC())
        pbc_odd = LineLattice(5, PBC())
        obc = LineLattice(5, OBC())

        @test topology(pbc_even) isa TopologyTrait{:line}
        @test periodicity(pbc_even) isa Periodic
        @test periodicity(obc) isa Aperiodic

        @test is_bipartite(pbc_even) == true    # even cycle
        @test is_bipartite(pbc_odd) == false    # odd cycle
        @test is_bipartite(obc) == true         # chain is always bipartite

        @test reciprocal_support(pbc_even) isa HasReciprocal
        @test reciprocal_support(obc) isa NoReciprocal
    end

    @testset "degenerate N=2 PBC (duplicate neighbour collapse)" begin
        lat = LineLattice(2, PBC())
        # Both left and right wrap to the other site → just one neighbour.
        @test neighbors(lat, 1) == [2]
        @test neighbors(lat, 2) == [1]
        @test length(collect(bonds(lat))) == 1
    end
end
