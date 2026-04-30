using LatticeCore
using Test

# Round-trip helper: site_index ∘ lattice_coord = identity on
# 1:prod(dims) * nsub for the given indexing.
function _assert_roundtrip_3d(indexing, dims, nsub)
    total = prod(dims) * nsub
    for i in 1:total
        lc = lattice_coord(indexing, dims, nsub, i)
        @test site_index(indexing, dims, nsub, lc) == i
    end
end

# Unique-cover helper: every (cell, sub) maps to a distinct index in
# 1:prod(dims) * nsub.
function _assert_unique_cover_3d(indexing, dims, nsub)
    seen = Set{Int}()
    for cx in 1:dims[1], cy in 1:dims[2], cz in 1:dims[3], s in 1:nsub
        push!(
            seen, site_index(indexing, dims, nsub, LatticeCoord((cx, cy, cz), s))
        )
    end
    @test length(seen) == prod(dims) * nsub
    @test extrema(seen) == (1, prod(dims) * nsub)
end

@testset "3D ColMajor indexing" begin
    @testset "site_index 2x2x2" begin
        # Convention: x slowest, z fastest in storage order.
        #   (1,1,1)->1, (1,1,2)->2, (1,2,1)->3, (1,2,2)->4,
        #   (2,1,1)->5, (2,1,2)->6, (2,2,1)->7, (2,2,2)->8.
        @test site_index(ColMajor(), (2, 2, 2), 1, LatticeCoord((1, 1, 1))) == 1
        @test site_index(ColMajor(), (2, 2, 2), 1, LatticeCoord((1, 1, 2))) == 2
        @test site_index(ColMajor(), (2, 2, 2), 1, LatticeCoord((1, 2, 1))) == 3
        @test site_index(ColMajor(), (2, 2, 2), 1, LatticeCoord((1, 2, 2))) == 4
        @test site_index(ColMajor(), (2, 2, 2), 1, LatticeCoord((2, 1, 1))) == 5
        @test site_index(ColMajor(), (2, 2, 2), 1, LatticeCoord((2, 2, 2))) == 8
    end

    @testset "lattice_coord 2x2x2" begin
        @test lattice_coord(ColMajor(), (2, 2, 2), 1, 1).cell == (1, 1, 1)
        @test lattice_coord(ColMajor(), (2, 2, 2), 1, 2).cell == (1, 1, 2)
        @test lattice_coord(ColMajor(), (2, 2, 2), 1, 3).cell == (1, 2, 1)
        @test lattice_coord(ColMajor(), (2, 2, 2), 1, 5).cell == (2, 1, 1)
        @test lattice_coord(ColMajor(), (2, 2, 2), 1, 8).cell == (2, 2, 2)
    end

    @testset "round-trip 3D ColMajor" begin
        _assert_roundtrip_3d(ColMajor(), (2, 2, 2), 1)
        _assert_roundtrip_3d(ColMajor(), (3, 3, 3), 1)
        _assert_roundtrip_3d(ColMajor(), (3, 2, 4), 1)
        _assert_roundtrip_3d(ColMajor(), (3, 3, 3), 2)
    end

    @testset "every site index is hit exactly once (3D ColMajor)" begin
        for dims in ((2, 2, 2), (3, 3, 3), (3, 2, 4)), nsub in (1, 2)
            _assert_unique_cover_3d(ColMajor(), dims, nsub)
        end
    end

    @testset "sublattice is innermost" begin
        # nsub = 2: index advances by sublattice first, then by cz.
        @test site_index(ColMajor(), (2, 2, 2), 2, LatticeCoord((1, 1, 1), 1)) == 1
        @test site_index(ColMajor(), (2, 2, 2), 2, LatticeCoord((1, 1, 1), 2)) == 2
        @test site_index(ColMajor(), (2, 2, 2), 2, LatticeCoord((1, 1, 2), 1)) == 3
    end
end

@testset "3D Snake indexing" begin
    @testset "site_index 3x3x3 boustrophedon" begin
        # Layer z=1 (forward y): rows y=1 forward, y=2 reverse, y=3 forward.
        #   (1,1,1)->1, (2,1,1)->2, (3,1,1)->3,
        #   (3,2,1)->4, (2,2,1)->5, (1,2,1)->6,
        #   (1,3,1)->7, (2,3,1)->8, (3,3,1)->9.
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 1, 1))) == 1
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 1, 1))) == 3
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 2, 1))) == 4
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 2, 1))) == 6
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 3, 1))) == 7
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 3, 1))) == 9

        # Layer z=2 (reverse y): rows traverse y=3, y=2, y=1.
        # effective_y for cy=3 is 1 (odd → x forward), so (1,3,2)->10.
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 3, 2))) == 10
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 3, 2))) == 12
        # effective_y for cy=2 is 2 (even → x reversed), so (3,2,2)->13.
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 2, 2))) == 13
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 2, 2))) == 15
        # effective_y for cy=1 is 3 (odd → x forward), so (1,1,2)->16.
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 1, 2))) == 16
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 1, 2))) == 18

        # Layer z=3 (forward y again): start back at (1,1,3).
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((1, 1, 3))) == 19
        @test site_index(Snake(), (3, 3, 3), 1, LatticeCoord((3, 3, 3))) == 27
    end

    @testset "lattice_coord 3x3x3" begin
        # Inverse round-trip on a few representative indices.
        @test lattice_coord(Snake(), (3, 3, 3), 1, 1).cell == (1, 1, 1)
        @test lattice_coord(Snake(), (3, 3, 3), 1, 4).cell == (3, 2, 1)
        @test lattice_coord(Snake(), (3, 3, 3), 1, 9).cell == (3, 3, 1)
        @test lattice_coord(Snake(), (3, 3, 3), 1, 10).cell == (1, 3, 2)
        @test lattice_coord(Snake(), (3, 3, 3), 1, 18).cell == (3, 1, 2)
        @test lattice_coord(Snake(), (3, 3, 3), 1, 19).cell == (1, 1, 3)
        @test lattice_coord(Snake(), (3, 3, 3), 1, 27).cell == (3, 3, 3)
    end

    @testset "round-trip 3D Snake" begin
        _assert_roundtrip_3d(Snake(), (3, 3, 3), 1)
        _assert_roundtrip_3d(Snake(), (2, 2, 2), 1)
        _assert_roundtrip_3d(Snake(), (4, 3, 2), 1)
        _assert_roundtrip_3d(Snake(), (3, 3, 3), 2)
        _assert_roundtrip_3d(Snake(), (4, 3, 2), 3)
    end

    @testset "every site index is hit exactly once (3D Snake)" begin
        for dims in ((2, 2, 2), (3, 3, 3), (4, 3, 2)), nsub in (1, 2)
            _assert_unique_cover_3d(Snake(), dims, nsub)
        end
    end

    @testset "consecutive indices are physically adjacent" begin
        # The defining property of a snake/boustrophedon ordering: any
        # two consecutive site indices live on (Manhattan-)adjacent
        # cells. This catches off-by-one bugs in the parity/reverse
        # logic at layer transitions.
        for dims in ((3, 3, 3), (4, 3, 2), (2, 4, 3))
            total = prod(dims)
            for i in 1:(total - 1)
                lc1 = lattice_coord(Snake(), dims, 1, i)
                lc2 = lattice_coord(Snake(), dims, 1, i + 1)
                d = sum(abs.(lc1.cell .- lc2.cell))
                @test d == 1
            end
        end
    end
end

@testset "3D RowMajor / ColMajor / Snake produce distinct permutations" begin
    # All three indexings must cover the same set of sites (so the
    # set of indices is identical), but the orderings must differ.
    dims = (3, 3, 3)
    nsub = 2
    total = prod(dims) * nsub

    function _coord_sequence(indexing)
        return [lattice_coord(indexing, dims, nsub, i) for i in 1:total]
    end

    row = _coord_sequence(RowMajor())
    col = _coord_sequence(ColMajor())
    snake = _coord_sequence(Snake())

    # Distinct permutations.
    @test row != col
    @test row != snake
    @test col != snake

    # Same underlying set.
    @test Set(row) == Set(col) == Set(snake)
    @test length(Set(row)) == total
end
