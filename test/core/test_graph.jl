using LatticeCore
using LinearAlgebra
using SparseArrays
using StaticArrays
using Test

# Minimal local lattice with custom topology for testing graph ops
# without depending on the reference lattices' particular bond order.

# A path graph 1-2-3-4-5 (1D chain, OBC-like).
struct GraphPath <: AbstractLattice{1,Float64}
    N::Int
end
LatticeCore.num_sites(l::GraphPath) = l.N
LatticeCore.position(l::GraphPath, i::Int) = SVector{1,Float64}(Float64(i))
function LatticeCore.neighbors(l::GraphPath, i::Int)
    return filter(j -> 1 <= j <= l.N, [i - 1, i + 1])
end
LatticeCore.boundary(::GraphPath) = nothing
LatticeCore.size_trait(l::GraphPath) = FiniteSize((l.N,))

# Two disjoint triangles: {1,2,3} and {4,5,6}.
struct TwoTriangles <: AbstractLattice{1,Float64} end
LatticeCore.num_sites(::TwoTriangles) = 6
LatticeCore.position(::TwoTriangles, i::Int) = SVector{1,Float64}(Float64(i))
function LatticeCore.neighbors(::TwoTriangles, i::Int)
    return if i == 1
        [2, 3]
    elseif i == 2
        [1, 3]
    elseif i == 3
        [1, 2]
    elseif i == 4
        [5, 6]
    elseif i == 5
        [4, 6]
    elseif i == 6
        [4, 5]
    else
        Int[]
    end
end
LatticeCore.boundary(::TwoTriangles) = nothing
LatticeCore.size_trait(::TwoTriangles) = FiniteSize((6,))

@testset "Graph operations" begin
    @testset "adjacency_matrix on path graph" begin
        lat = GraphPath(5)
        A = adjacency_matrix(lat)
        @test A isa SparseMatrixCSC{Bool,Int}
        @test size(A) == (5, 5)
        @test A == transpose(A)
        # 4 undirected edges, 8 stored true entries
        @test count(A) == 8
        @test A[1, 2] && A[2, 1]
        @test A[2, 3] && A[3, 4] && A[4, 5]
        @test !A[1, 3]
        @test !A[1, 1]   # no self-loops

        Ad = adjacency_matrix(lat; sparse=false)
        @test Ad isa Matrix{Bool}
        @test Ad == Matrix(A)
    end

    @testset "adjacency_matrix on SimpleSquareLattice" begin
        sq = SimpleSquareLattice(3, 3, OpenAxis())
        A = adjacency_matrix(sq)
        @test size(A) == (9, 9)
        @test A == transpose(A)
        # 3x3 OBC square has 12 edges -> 24 stored true entries.
        @test count(A) == 24
    end

    @testset "shortest_path on path graph" begin
        lat = GraphPath(5)
        d, p = shortest_path(lat, 1, 5)
        @test d == 4
        @test p == [1, 2, 3, 4, 5]

        d2, p2 = shortest_path(lat, 3, 3)
        @test d2 == 0
        @test p2 == [3]

        d3, p3 = shortest_path(lat, 4, 2)
        @test d3 == 2
        @test p3 == [4, 3, 2]
    end

    @testset "shortest_path on disconnected graph" begin
        lat = TwoTriangles()
        d, p = shortest_path(lat, 1, 4)
        @test d == typemax(Int)
        @test isempty(p)

        # Within one component
        d2, p2 = shortest_path(lat, 4, 6)
        @test d2 == 1
        @test p2 == [4, 6]
    end

    @testset "shortest_path bounds checking" begin
        lat = GraphPath(3)
        @test_throws BoundsError shortest_path(lat, 0, 2)
        @test_throws BoundsError shortest_path(lat, 1, 4)
    end

    @testset "shortest_path on SimpleSquareLattice OBC" begin
        sq = SimpleSquareLattice(3, 3, OpenAxis())
        d, p = shortest_path(sq, 1, 9)
        @test d == 4
        @test first(p) == 1 && last(p) == 9
        @test length(p) == 5
        # All consecutive pairs in the path must be neighbours.
        for k in 1:(length(p) - 1)
            @test p[k + 1] in neighbors(sq, p[k])
        end
    end

    @testset "connected_components" begin
        lat1 = GraphPath(5)
        comps1 = connected_components(lat1)
        @test length(comps1) == 1
        @test comps1[1] == [1, 2, 3, 4, 5]

        lat2 = TwoTriangles()
        comps2 = connected_components(lat2)
        @test length(comps2) == 2
        @test comps2[1] == [1, 2, 3]
        @test comps2[2] == [4, 5, 6]
    end

    @testset "shortest_path with uniform weights matches BFS" begin
        # When `weights === identity_weight` (the default), the API
        # falls back to BFS and returns the Int hop count — see the
        # `shortest_path` docstring. Passing `identity_weight` by
        # keyword should be observationally identical to the no-kwarg
        # call shape.
        for lat in (GraphPath(6), SimpleSquareLattice(3, 3, OpenAxis()))
            N = num_sites(lat)
            for s in 1:N, t in 1:N
                d_bfs, p_bfs = shortest_path(lat, s, t)
                d_kw, p_kw = shortest_path(lat, s, t; weights=identity_weight)
                @test d_kw === d_bfs
                @test p_kw == p_bfs
            end
        end

        # A user-supplied uniform-1 callback (≠ identity_weight by
        # `===`) takes the Dijkstra path and must produce the same
        # connectivity-distance up to Float64 promotion.
        uniform_w(_lat, _b) = 1.0
        for lat in (GraphPath(6), SimpleSquareLattice(3, 3, OpenAxis()))
            N = num_sites(lat)
            for s in 1:N, t in 1:N
                d_bfs, p_bfs = shortest_path(lat, s, t)
                c_dij, p_dij = shortest_path(lat, s, t; weights=uniform_w)
                if d_bfs == typemax(Int)
                    @test c_dij == Inf
                    @test isempty(p_dij)
                else
                    @test c_dij == Float64(d_bfs)
                    @test first(p_dij) == s && last(p_dij) == t
                    @test length(p_dij) == length(p_bfs)
                end
            end
        end
    end

    @testset "shortest_path with custom weights" begin
        # Path graph 1-2-3-4-5; bond (i, i+1) weight = 1 if i odd, 2 if i even.
        lat = GraphPath(5)
        alt_w(_lat, b) = isodd(min(b.i, b.j)) ? 1.0 : 2.0

        c, p = shortest_path(lat, 1, 5; weights=alt_w)
        # Costs along 1-2-3-4-5: 1 + 2 + 1 + 2 = 6
        @test c == 6.0
        @test p == [1, 2, 3, 4, 5]

        # 2 -> 5 traverses bonds (2,3), (3,4), (4,5) with weights 2, 1, 2.
        c2, p2 = shortest_path(lat, 2, 5; weights=alt_w)
        @test c2 == 5.0
        @test p2 == [2, 3, 4, 5]

        # src == dst stays zero-cost regardless of weights.
        c0, p0 = shortest_path(lat, 3, 3; weights=alt_w)
        @test c0 == 0.0
        @test p0 == [3]
    end

    @testset "shortest_path on disconnected graph (weighted)" begin
        lat = TwoTriangles()
        # A genuine custom callback engages the Dijkstra path; cross-
        # component queries should return `(Inf, Int[])`.
        uniform_w(_lat, _b) = 1.0
        c, p = shortest_path(lat, 1, 4; weights=uniform_w)
        @test c == Inf
        @test isempty(p)
        # Heavy custom weights stay Inf across the cut too.
        big_w(_lat, _b) = 100.0
        c2, p2 = shortest_path(lat, 2, 5; weights=big_w)
        @test c2 == Inf
        @test isempty(p2)
        # Within a component, total cost = 100 * (one hop) = 100.
        c3, p3 = shortest_path(lat, 4, 6; weights=big_w)
        @test c3 == 100.0
        @test first(p3) == 4 && last(p3) == 6
    end

    @testset "shortest_path keyword bounds checking" begin
        lat = GraphPath(3)
        @test_throws BoundsError shortest_path(lat, 0, 2; weights=identity_weight)
        @test_throws BoundsError shortest_path(lat, 1, 4; weights=identity_weight)
    end

    @testset "distance_matrix on small SimpleSquareLattice" begin
        sq = SimpleSquareLattice(4, 4, OpenAxis())
        D = distance_matrix(sq)
        @test D isa Matrix{Float64}
        @test size(D) == (16, 16)
        @test D == transpose(D)              # undirected -> symmetric
        @test all(iszero, diag(D))
        # Spot-check: corner-to-corner Manhattan distance on a 4x4
        # OBC square is 6 hops.
        @test D[1, 16] == 6.0
        # Same answer via shortest_path
        d_sp, _ = shortest_path(sq, 1, 16; weights=identity_weight)
        @test d_sp == D[1, 16]
    end

    @testset "distance_matrix Floyd-Warshall vs Dijkstra agree" begin
        # Force the threshold to flip algorithms on the same lattice
        # and check both code paths produce the same matrix.
        lat = GraphPath(8)
        D_fw = distance_matrix(lat; floyd_warshall_threshold=100)
        D_dj = distance_matrix(lat; floyd_warshall_threshold=0)
        @test D_fw == D_dj
        # And the matrix entries must agree element-wise with the
        # single-pair shortest_path query (Float64-promoted to match
        # `distance_matrix`'s element type).
        for s in 1:num_sites(lat), t in 1:num_sites(lat)
            d_sp, _ = shortest_path(lat, s, t)
            expected = d_sp == typemax(Int) ? Inf : Float64(d_sp)
            @test D_fw[s, t] == expected
        end
    end

    @testset "distance_matrix on disconnected graph" begin
        lat = TwoTriangles()
        D = distance_matrix(lat)
        @test size(D) == (6, 6)
        @test all(iszero, diag(D))
        # Cross-component entries are Inf
        @test D[1, 4] == Inf
        @test D[3, 5] == Inf
        # Within-component entries are finite
        @test D[1, 2] == 1.0
        @test D[4, 6] == 1.0
        # Symmetry
        @test D == transpose(D)
    end

    @testset "distance_matrix with custom weights" begin
        lat = GraphPath(4)
        # Bond (i, i+1) weight = i (so 1, 2, 3 along the chain).
        ramp_w(_lat, b) = Float64(min(b.i, b.j))
        D = distance_matrix(lat; weights=ramp_w)
        # 1->2: 1, 1->3: 1+2=3, 1->4: 1+2+3=6
        @test D[1, 2] == 1.0
        @test D[1, 3] == 3.0
        @test D[1, 4] == 6.0
        # Symmetric
        @test D == transpose(D)
    end
end
