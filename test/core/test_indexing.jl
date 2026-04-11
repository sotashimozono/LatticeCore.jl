using LatticeCore
using Test

# Exhaustive round-trip helper: checks site_index ∘ lattice_coord = identity
# on every valid site index for a given indexing / dims / nsub.
function _assert_roundtrip(indexing, dims, nsub)
    total = prod(dims) * nsub
    for i in 1:total
        lc = lattice_coord(indexing, dims, nsub, i)
        @test site_index(indexing, dims, nsub, lc) == i
    end
end

@testset "AbstractIndexing" begin
    @testset "type hierarchy" begin
        @test RowMajor <: AbstractIndexing
        @test ColMajor <: AbstractIndexing
        @test Snake <: AbstractIndexing
    end

    @testset "RowMajor 1D" begin
        # nsub = 1
        @test site_index(RowMajor(), (5,), 1, LatticeCoord((1,))) == 1
        @test site_index(RowMajor(), (5,), 1, LatticeCoord((5,))) == 5
        @test lattice_coord(RowMajor(), (5,), 1, 3).cell == (3,)

        _assert_roundtrip(RowMajor(), (5,), 1)

        # nsub = 2 (two sublattices per cell)
        @test site_index(RowMajor(), (3,), 2, LatticeCoord((1,), 1)) == 1
        @test site_index(RowMajor(), (3,), 2, LatticeCoord((1,), 2)) == 2
        @test site_index(RowMajor(), (3,), 2, LatticeCoord((2,), 1)) == 3
        @test site_index(RowMajor(), (3,), 2, LatticeCoord((3,), 2)) == 6

        _assert_roundtrip(RowMajor(), (3,), 2)
    end

    @testset "RowMajor 2D" begin
        # 3x2 row-major layout:
        #   row 1 (y=1): sites 1, 2, 3
        #   row 2 (y=2): sites 4, 5, 6
        @test site_index(RowMajor(), (3, 2), 1, LatticeCoord((1, 1))) == 1
        @test site_index(RowMajor(), (3, 2), 1, LatticeCoord((3, 1))) == 3
        @test site_index(RowMajor(), (3, 2), 1, LatticeCoord((1, 2))) == 4
        @test site_index(RowMajor(), (3, 2), 1, LatticeCoord((3, 2))) == 6

        @test lattice_coord(RowMajor(), (3, 2), 1, 4).cell == (1, 2)
        @test lattice_coord(RowMajor(), (3, 2), 1, 6).cell == (3, 2)

        _assert_roundtrip(RowMajor(), (3, 2), 1)
        _assert_roundtrip(RowMajor(), (4, 5), 1)
        _assert_roundtrip(RowMajor(), (3, 2), 2)   # two sublattices
    end

    @testset "ColMajor 2D" begin
        # 3x2 col-major layout:
        #   col 1 (x=1): sites 1, 2
        #   col 2 (x=2): sites 3, 4
        #   col 3 (x=3): sites 5, 6
        @test site_index(ColMajor(), (3, 2), 1, LatticeCoord((1, 1))) == 1
        @test site_index(ColMajor(), (3, 2), 1, LatticeCoord((1, 2))) == 2
        @test site_index(ColMajor(), (3, 2), 1, LatticeCoord((2, 1))) == 3
        @test site_index(ColMajor(), (3, 2), 1, LatticeCoord((3, 2))) == 6

        _assert_roundtrip(ColMajor(), (3, 2), 1)
        _assert_roundtrip(ColMajor(), (4, 5), 1)
        _assert_roundtrip(ColMajor(), (3, 2), 3)
    end

    @testset "Snake 2D" begin
        # 3x3 snake layout:
        #   y=1 (forward): 1, 2, 3
        #   y=2 (reverse): 4, 5, 6  ↔  (3,2), (2,2), (1,2)
        #   y=3 (forward): 7, 8, 9
        @test site_index(Snake(), (3, 3), 1, LatticeCoord((1, 1))) == 1
        @test site_index(Snake(), (3, 3), 1, LatticeCoord((3, 1))) == 3
        @test site_index(Snake(), (3, 3), 1, LatticeCoord((3, 2))) == 4
        @test site_index(Snake(), (3, 3), 1, LatticeCoord((1, 2))) == 6
        @test site_index(Snake(), (3, 3), 1, LatticeCoord((1, 3))) == 7
        @test site_index(Snake(), (3, 3), 1, LatticeCoord((3, 3))) == 9

        # Inverse
        @test lattice_coord(Snake(), (3, 3), 1, 4).cell == (3, 2)
        @test lattice_coord(Snake(), (3, 3), 1, 6).cell == (1, 2)

        _assert_roundtrip(Snake(), (3, 3), 1)
        _assert_roundtrip(Snake(), (4, 3), 1)
        _assert_roundtrip(Snake(), (3, 4), 2)
    end

    @testset "round-trip count covers every site" begin
        # Ensure every site index is the image of exactly one LatticeCoord
        # for each indexing. This rules out collisions and holes.
        for indexing in (RowMajor(), ColMajor(), Snake())
            for dims in ((3, 3), (4, 5), (5, 2))
                for nsub in (1, 2)
                    seen = Set{Int}()
                    for cx in 1:dims[1], cy in 1:dims[2], s in 1:nsub
                        push!(
                            seen,
                            site_index(indexing, dims, nsub, LatticeCoord((cx, cy), s)),
                        )
                    end
                    @test length(seen) == prod(dims) * nsub
                    @test extrema(seen) == (1, prod(dims) * nsub)
                end
            end
        end
    end
end
