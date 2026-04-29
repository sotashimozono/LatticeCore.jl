using LatticeCore
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
end
