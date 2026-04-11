using LatticeCore
using StaticArrays
using Test

@testset "LineLattice" begin
    @testset "construction (default Float64, default PBC)" begin
        lat = LineLattice(5)
        @test lat isa LineLattice{Float64}
        @test lat isa AbstractLattice{1,Float64}
        @test dimension(lat) == 1
        @test eltype(lat) === Float64
        @test num_sites(lat) == 5

        bc = boundary(lat)
        @test bc isa LatticeBoundary{1}
        @test bc.axes[1] isa PeriodicAxis
        @test bc.modifier isa NoModifier
    end

    @testset "axis BC shorthand construction" begin
        @test boundary(LineLattice(5, PeriodicAxis())).axes[1] isa PeriodicAxis
        @test boundary(LineLattice(5, OpenAxis())).axes[1] isa OpenAxis
        @test boundary(LineLattice(5, TwistedAxis(0.3))).axes[1] isa TwistedAxis
    end

    @testset "explicit LatticeBoundary construction" begin
        lat = LineLattice(5, LatticeBoundary((PeriodicAxis(),), NoModifier()))
        @test boundary(lat).axes[1] isa PeriodicAxis
    end

    @testset "PeriodicAxis neighbors" begin
        lat = LineLattice(5, PeriodicAxis())
        @test Set(neighbors(lat, 3)) == Set([2, 4])
        @test Set(neighbors(lat, 1)) == Set([5, 2])
        @test Set(neighbors(lat, 5)) == Set([4, 1])
    end

    @testset "OpenAxis neighbors" begin
        lat = LineLattice(5, OpenAxis())
        @test neighbors(lat, 1) == [2]
        @test Set(neighbors(lat, 3)) == Set([2, 4])
        @test neighbors(lat, 5) == [4]
    end

    @testset "positions" begin
        lat = LineLattice(4, OpenAxis())
        ps = collect(positions(lat))
        @test length(ps) == 4
        @test ps[1] == SVector(1.0)
        @test ps[4] == SVector(4.0)
    end

    @testset "bond counts" begin
        # OBC chain of length N has N-1 bonds
        @test length(collect(bonds(LineLattice(5, OpenAxis())))) == 4
        # PBC cycle of length N has N bonds
        @test length(collect(bonds(LineLattice(5, PeriodicAxis())))) == 5
    end

    @testset "neighbor_bonds with wrapped unit displacements" begin
        lat = LineLattice(5, PeriodicAxis())

        # Interior site: both bonds have ±1 displacement.
        nb_mid = collect(neighbor_bonds(lat, 3))
        @test length(nb_mid) == 2
        @test all(b.i == 3 for b in nb_mid)
        @test Set(b.j for b in nb_mid) == Set([2, 4])
        @test Set(b.vector[1] for b in nb_mid) == Set([1.0, -1.0])

        # Boundary-crossing site: +1 direction wraps to site 1 (not raw +4).
        nb_edge = collect(neighbor_bonds(lat, 5))
        @test length(nb_edge) == 2
        @test Set(b.j for b in nb_edge) == Set([4, 1])
        # The bond (5 -> 1) must carry a +1 displacement, not a raw +(-4).
        wrap_bond = first(b for b in nb_edge if b.j == 1)
        @test wrap_bond.vector == SVector(1.0)
    end

    @testset "size trait / is_finite" begin
        lat = LineLattice(7, OpenAxis())
        @test size_trait(lat) isa FiniteSize{1}
        @test size_trait(lat).dims == (7,)
        @test is_finite(lat) == true
    end

    @testset "trait overrides" begin
        pbc_even = LineLattice(6, PeriodicAxis())
        pbc_odd = LineLattice(5, PeriodicAxis())
        obc = LineLattice(5, OpenAxis())
        twisted = LineLattice(6, TwistedAxis(0.2))

        @test topology(pbc_even) isa TopologyTrait{:line}
        @test periodicity(pbc_even) isa Periodic
        @test periodicity(obc) isa Aperiodic
        @test periodicity(twisted) isa Periodic     # twisted ⇒ still cyclic

        @test is_bipartite(pbc_even) == true       # even cycle
        @test is_bipartite(pbc_odd) == false       # odd cycle
        @test is_bipartite(obc) == true            # chain is always bipartite
        @test is_bipartite(twisted) == true        # even twisted cycle

        @test reciprocal_support(pbc_even) isa HasReciprocal
        @test reciprocal_support(obc) isa NoReciprocal
        @test reciprocal_support(twisted) isa HasReciprocal
    end

    @testset "degenerate N=2 PBC (duplicate neighbour collapse)" begin
        lat = LineLattice(2, PeriodicAxis())
        @test neighbors(lat, 1) == [2]
        @test neighbors(lat, 2) == [1]
        @test length(collect(bonds(lat))) == 1
    end

    @testset "coordinate conversions" begin
        lat = LineLattice(5, PeriodicAxis())

        rs = to_real(lat, LatticeCoord((3,)))
        @test rs isa RealSpace{1,Float64}
        @test rs.x == SVector(3.0)

        lc = to_lattice(lat, RealSpace((4.0,)))
        @test lc isa LatticeCoord{1}
        @test lc.cell == (4,)
        @test lc.sublattice == 1

        # Round-trip
        for i in 1:num_sites(lat)
            lc_i = LatticeCoord((i,))
            @test to_lattice(lat, to_real(lat, lc_i)) == lc_i
        end
    end
end
