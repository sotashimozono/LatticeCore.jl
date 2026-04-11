using LatticeCore
using StaticArrays
using Test

@testset "SimpleSquareLattice" begin
    @testset "construction" begin
        lat = SimpleSquareLattice(3, 4)
        @test lat isa SimpleSquareLattice{Float64}
        @test lat isa AbstractLattice{2,Float64}
        @test dimension(lat) == 2
        @test num_sites(lat) == 12

        bc = boundary(lat)
        @test bc isa LatticeBoundary{2}
        @test bc.axes[1] isa PeriodicAxis
        @test bc.axes[2] isa PeriodicAxis
    end

    @testset "axis BC shorthand" begin
        sq = SimpleSquareLattice(3, 3, OpenAxis())
        @test boundary(sq).axes[1] isa OpenAxis
        @test boundary(sq).axes[2] isa OpenAxis
    end

    @testset "explicit mixed BC (cylinder)" begin
        cyl = SimpleSquareLattice(3, 3, LatticeBoundary((PeriodicAxis(), OpenAxis())))
        @test boundary(cyl).axes[1] isa PeriodicAxis
        @test boundary(cyl).axes[2] isa OpenAxis
    end

    @testset "position layout (row-major)" begin
        lat = SimpleSquareLattice(3, 2, OpenAxis())
        @test position(lat, 1) == SVector(1.0, 1.0)
        @test position(lat, 2) == SVector(2.0, 1.0)
        @test position(lat, 3) == SVector(3.0, 1.0)
        @test position(lat, 4) == SVector(1.0, 2.0)
        @test position(lat, 5) == SVector(2.0, 2.0)
        @test position(lat, 6) == SVector(3.0, 2.0)
    end

    @testset "PBC neighbors" begin
        lat = SimpleSquareLattice(3, 3, PeriodicAxis())
        # Interior site (2,2) = site 5
        @test Set(neighbors(lat, 5)) == Set([6, 4, 8, 2])
        # Corner site (1,1) = site 1; wraps to (3,1) and (1,3)
        @test Set(neighbors(lat, 1)) == Set([2, 3, 4, 7])
    end

    @testset "OBC neighbors" begin
        lat = SimpleSquareLattice(3, 3, OpenAxis())
        @test Set(neighbors(lat, 5)) == Set([6, 4, 8, 2])
        @test Set(neighbors(lat, 1)) == Set([2, 4])
        @test Set(neighbors(lat, 2)) == Set([1, 3, 5])
        @test Set(neighbors(lat, 9)) == Set([8, 6])
    end

    @testset "cylinder neighbors (x PBC, y OBC)" begin
        cyl = SimpleSquareLattice(3, 3, LatticeBoundary((PeriodicAxis(), OpenAxis())))
        # Corner (1,1) = site 1: x wraps to (3,1), y open so only (1,2)
        @test Set(neighbors(cyl, 1)) == Set([2, 3, 4])
        # Middle (2,2) = site 5: full 4 neighbours
        @test Set(neighbors(cyl, 5)) == Set([6, 4, 8, 2])
        # Top row (1,3) = site 7: wraps in x, open in y → (2,3), (3,3), (1,2)
        @test Set(neighbors(cyl, 7)) == Set([8, 9, 4])
    end

    @testset "bond counts" begin
        # PBC Lx x Ly: 2 * Lx * Ly bonds
        @test length(collect(bonds(SimpleSquareLattice(3, 3, PeriodicAxis())))) == 18
        @test length(collect(bonds(SimpleSquareLattice(4, 5, PeriodicAxis())))) == 40

        # OBC Lx x Ly: (Lx-1) * Ly + Lx * (Ly-1)
        @test length(collect(bonds(SimpleSquareLattice(3, 3, OpenAxis())))) == 12
        @test length(collect(bonds(SimpleSquareLattice(4, 5, OpenAxis())))) ==
            (3 * 5) + (4 * 4)

        # Cylinder: Lx PBC bonds per row (Lx) * Ly rows  +  (Ly-1) vertical bonds * Lx columns
        cyl = SimpleSquareLattice(3, 3, LatticeBoundary((PeriodicAxis(), OpenAxis())))
        @test length(collect(bonds(cyl))) == (3 * 3) + (2 * 3)   # 9 + 6 = 15
    end

    @testset "bond wrapping: PBC corner bond carries unit vector" begin
        lat = SimpleSquareLattice(3, 3, PeriodicAxis())
        # Site 1 = (1,1): -x neighbour wraps to site 3 = (3,1).
        # The bond vector must be (-1, 0), not the raw (2, 0).
        nb1 = collect(neighbor_bonds(lat, 1))
        wrap_x = first(b for b in nb1 if b.j == 3)
        @test wrap_x.vector == SVector(-1.0, 0.0)

        # Similarly for the y-wrapping bond to site 7 = (1, 3).
        wrap_y = first(b for b in nb1 if b.j == 7)
        @test wrap_y.vector == SVector(0.0, -1.0)
    end

    @testset "size trait / is_finite" begin
        lat = SimpleSquareLattice(4, 6, PeriodicAxis())
        st = size_trait(lat)
        @test st isa FiniteSize{2}
        @test st.dims == (4, 6)
        @test is_finite(lat) == true
    end

    @testset "trait overrides" begin
        pbc_even = SimpleSquareLattice(4, 4, PeriodicAxis())
        pbc_odd = SimpleSquareLattice(3, 3, PeriodicAxis())
        obc = SimpleSquareLattice(3, 3, OpenAxis())
        cyl = SimpleSquareLattice(4, 3, LatticeBoundary((PeriodicAxis(), OpenAxis())))

        @test topology(pbc_even) isa TopologyTrait{:square}
        @test periodicity(pbc_even) isa Periodic
        @test periodicity(obc) isa Aperiodic
        @test periodicity(cyl) isa Aperiodic   # cylinder has an open axis

        # Bipartiteness
        @test is_bipartite(pbc_even) == true
        @test is_bipartite(pbc_odd) == false
        @test is_bipartite(SimpleSquareLattice(4, 3, PeriodicAxis())) == false
        @test is_bipartite(obc) == true
        @test is_bipartite(cyl) == true           # PBC axis is even (Lx=4)

        @test reciprocal_support(pbc_even) isa HasReciprocal
        @test reciprocal_support(obc) isa NoReciprocal
        @test reciprocal_support(cyl) isa NoReciprocal  # mixed ⇒ no uniform recip
    end

    @testset "bond_center via default Bond iteration" begin
        lat = SimpleSquareLattice(3, 3, OpenAxis())
        b = Bond{2,Float64}(1, 2, SVector(1.0, 0.0), :nearest)
        @test bond_center(lat, b) == SVector(1.5, 1.0)
    end
end
