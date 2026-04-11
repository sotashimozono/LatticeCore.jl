using LatticeCore
using StaticArrays
using Test

@testset "SimpleSquareLattice" begin
    @testset "construction" begin
        lat = SimpleSquareLattice(3, 4)
        @test lat isa SimpleSquareLattice{Float64,PBC}
        @test lat isa AbstractLattice{2,Float64}
        @test dimension(lat) == 2
        @test num_sites(lat) == 12
        @test boundary(lat) isa PBC
    end

    @testset "position layout (row-major)" begin
        lat = SimpleSquareLattice(3, 2, OBC())
        # Row-major: sites 1,2,3 = row y=1; sites 4,5,6 = row y=2.
        @test position(lat, 1) == SVector(1.0, 1.0)
        @test position(lat, 2) == SVector(2.0, 1.0)
        @test position(lat, 3) == SVector(3.0, 1.0)
        @test position(lat, 4) == SVector(1.0, 2.0)
        @test position(lat, 5) == SVector(2.0, 2.0)
        @test position(lat, 6) == SVector(3.0, 2.0)
    end

    @testset "PBC neighbors (interior / edges)" begin
        lat = SimpleSquareLattice(3, 3, PBC())

        # Interior site (2,2) = site 5
        @test Set(neighbors(lat, 5)) == Set([
            6,  # (3,2)
            4,  # (1,2)
            8,  # (2,3)
            2,  # (2,1)
        ])

        # Corner site (1,1) = site 1; wraps to (3,1) and (1,3)
        @test Set(neighbors(lat, 1)) == Set([
            2,  # (2,1)
            3,  # (3,1) via wrap
            4,  # (1,2)
            7,  # (1,3) via wrap
        ])
    end

    @testset "OBC neighbors (interior / corner / edge)" begin
        lat = SimpleSquareLattice(3, 3, OBC())

        # Interior site (2,2) = site 5
        @test Set(neighbors(lat, 5)) == Set([6, 4, 8, 2])

        # Corner site (1,1) = site 1: no wrap
        @test Set(neighbors(lat, 1)) == Set([2, 4])

        # Edge site (2,1) = site 2: has x-left, x-right, y-up only
        @test Set(neighbors(lat, 2)) == Set([1, 3, 5])

        # Opposite corner (3,3) = site 9
        @test Set(neighbors(lat, 9)) == Set([8, 6])
    end

    @testset "bond counts" begin
        # PBC Lx x Ly: 2 * Lx * Ly bonds (horizontal + vertical)
        lat_pbc = SimpleSquareLattice(3, 3, PBC())
        @test length(collect(bonds(lat_pbc))) == 18

        lat_pbc_4x5 = SimpleSquareLattice(4, 5, PBC())
        @test length(collect(bonds(lat_pbc_4x5))) == 40

        # OBC Lx x Ly: (Lx-1) * Ly  +  Lx * (Ly-1) bonds
        lat_obc = SimpleSquareLattice(3, 3, OBC())
        @test length(collect(bonds(lat_obc))) == 12

        lat_obc_4x5 = SimpleSquareLattice(4, 5, OBC())
        @test length(collect(bonds(lat_obc_4x5))) == (3 * 5) + (4 * 4)
    end

    @testset "size trait / is_finite" begin
        lat = SimpleSquareLattice(4, 6, PBC())
        st = size_trait(lat)
        @test st isa FiniteSize{2}
        @test st.dims == (4, 6)
        @test is_finite(lat) == true
    end

    @testset "trait overrides" begin
        pbc_even = SimpleSquareLattice(4, 4, PBC())
        pbc_odd = SimpleSquareLattice(3, 3, PBC())
        obc = SimpleSquareLattice(3, 3, OBC())

        @test topology(pbc_even) isa TopologyTrait{:square}
        @test periodicity(pbc_even) isa Periodic
        @test periodicity(obc) isa Aperiodic

        # PBC square is bipartite iff both Lx, Ly even
        @test is_bipartite(pbc_even) == true
        @test is_bipartite(pbc_odd) == false
        @test is_bipartite(SimpleSquareLattice(4, 3, PBC())) == false
        @test is_bipartite(obc) == true

        @test reciprocal_support(pbc_even) isa HasReciprocal
        @test reciprocal_support(obc) isa NoReciprocal
    end

    @testset "bond_center via default Bond iteration" begin
        lat = SimpleSquareLattice(3, 3, OBC())
        # Pick the bond between sites 1 (1,1) and 2 (2,1).
        b = Bond{2,Float64}(1, 2, SVector(1.0, 0.0), :nearest)
        @test bond_center(lat, b) == SVector(1.5, 1.0)
    end
end
