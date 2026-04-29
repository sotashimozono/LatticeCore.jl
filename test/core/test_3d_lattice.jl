using LatticeCore
using StaticArrays
using SparseArrays
using Test

# Minimal 3D dummy lattice: 2 x 2 x 2 grid (8 sites) with OBC-style
# nearest-neighbour links. The intent is to exercise the
# `AbstractLattice{3, T}` interface + the new 3D RowMajor indexing,
# not to provide a production cubic lattice. Real concrete cubics
# (Simple Cubic / FCC / BCC) live in a future Lattice3D.jl.

struct ToyCubicLattice <: AbstractLattice{3,Float64}
    Lx::Int
    Ly::Int
    Lz::Int
end

ToyCubicLattice() = ToyCubicLattice(2, 2, 2)

LatticeCore.num_sites(l::ToyCubicLattice) = l.Lx * l.Ly * l.Lz

@inline function _to_xyz(l::ToyCubicLattice, i::Int)
    Lx = l.Lx
    Ly = l.Ly
    c = i - 1
    cx = c % Lx + 1
    cyz = c ÷ Lx
    cy = cyz % Ly + 1
    cz = cyz ÷ Ly + 1
    return (cx, cy, cz)
end

@inline function _from_xyz(l::ToyCubicLattice, x::Int, y::Int, z::Int)
    return ((z - 1) * l.Ly + (y - 1)) * l.Lx + (x - 1) + 1
end

function LatticeCore.position(l::ToyCubicLattice, i::Int)
    x, y, z = _to_xyz(l, i)
    return SVector{3,Float64}(Float64(x), Float64(y), Float64(z))
end

function LatticeCore.neighbors(l::ToyCubicLattice, i::Int)
    x, y, z = _to_xyz(l, i)
    out = Int[]
    for (dx, dy, dz) in
        ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
        nx, ny, nz = x + dx, y + dy, z + dz
        (1 <= nx <= l.Lx && 1 <= ny <= l.Ly && 1 <= nz <= l.Lz) || continue
        push!(out, _from_xyz(l, nx, ny, nz))
    end
    return out
end

LatticeCore.boundary(::ToyCubicLattice) = nothing
LatticeCore.size_trait(l::ToyCubicLattice) = FiniteSize((l.Lx, l.Ly, l.Lz))

@testset "AbstractLattice{3} interface (ToyCubicLattice)" begin
    lat = ToyCubicLattice(2, 2, 2)

    @testset "type parameters and basics" begin
        @test lat isa AbstractLattice{3,Float64}
        @test dimension(lat) == 3
        @test eltype(lat) === Float64
        @test num_sites(lat) == 8
    end

    @testset "positions are SVector{3, Float64}" begin
        @test position(lat, 1) isa SVector{3,Float64}
        @test position(lat, 1) == SVector{3,Float64}(1.0, 1.0, 1.0)
        @test position(lat, 8) == SVector{3,Float64}(2.0, 2.0, 2.0)
        ps = collect(positions(lat))
        @test length(ps) == 8
        @test all(p isa SVector{3,Float64} for p in ps)
    end

    @testset "neighbors and connectivity" begin
        # Site (1,1,1) has 3 neighbours in a 2x2x2 OBC cube.
        @test length(neighbors(lat, 1)) == 3
        # Total edges in a 2x2x2 OBC cube: 12. Sum of degrees == 24.
        total_deg = sum(length(neighbors(lat, i)) for i in 1:num_sites(lat))
        @test total_deg == 24
    end

    @testset "default bond iterator produces SVector{3} displacements" begin
        bs = collect(bonds(lat))
        @test length(bs) == 12
        @test all(b isa Bond{3,Float64} for b in bs)
        @test all(b.vector isa SVector{3,Float64} for b in bs)
    end

    @testset "size_trait reports 3D dims" begin
        st = size_trait(lat)
        @test st isa FiniteSize{3}
        @test st.dims == (2, 2, 2)
    end

    @testset "graph API works in 3D" begin
        A = adjacency_matrix(lat)
        @test A isa SparseMatrixCSC{Bool,Int}
        @test size(A) == (8, 8)
        @test A == transpose(A)
        @test count(A) == 24

        d, p = shortest_path(lat, 1, 8)
        # Diagonal of the 2x2x2 cube: distance is 3.
        @test d == 3
        @test first(p) == 1 && last(p) == 8

        comps = connected_components(lat)
        @test length(comps) == 1
        @test comps[1] == collect(1:8)
    end
end

@testset "3D RowMajor indexing" begin
    function _assert_roundtrip_3d(dims, nsub)
        total = prod(dims) * nsub
        for i in 1:total
            lc = lattice_coord(RowMajor(), dims, nsub, i)
            @test site_index(RowMajor(), dims, nsub, lc) == i
        end
    end

    @testset "site_index 2x2x2" begin
        # x is fastest, z is slowest:
        #   (1,1,1)->1, (2,1,1)->2, (1,2,1)->3, (2,2,1)->4,
        #   (1,1,2)->5, (2,1,2)->6, (1,2,2)->7, (2,2,2)->8.
        @test site_index(RowMajor(), (2, 2, 2), 1, LatticeCoord((1, 1, 1))) == 1
        @test site_index(RowMajor(), (2, 2, 2), 1, LatticeCoord((2, 1, 1))) == 2
        @test site_index(RowMajor(), (2, 2, 2), 1, LatticeCoord((1, 2, 1))) == 3
        @test site_index(RowMajor(), (2, 2, 2), 1, LatticeCoord((2, 2, 2))) == 8
    end

    @testset "lattice_coord 2x2x2" begin
        @test lattice_coord(RowMajor(), (2, 2, 2), 1, 1).cell == (1, 1, 1)
        @test lattice_coord(RowMajor(), (2, 2, 2), 1, 5).cell == (1, 1, 2)
        @test lattice_coord(RowMajor(), (2, 2, 2), 1, 8).cell == (2, 2, 2)
    end

    @testset "round-trip 3D" begin
        _assert_roundtrip_3d((2, 2, 2), 1)
        _assert_roundtrip_3d((3, 2, 4), 1)
        _assert_roundtrip_3d((2, 3, 2), 2)
    end

    @testset "every site index is hit exactly once (3D)" begin
        for dims in ((2, 2, 2), (3, 2, 4)), nsub in (1, 2)
            seen = Set{Int}()
            for cx in 1:dims[1], cy in 1:dims[2], cz in 1:dims[3], s in 1:nsub
                push!(seen, site_index(RowMajor(), dims, nsub, LatticeCoord((cx, cy, cz), s)))
            end
            @test length(seen) == prod(dims) * nsub
            @test extrema(seen) == (1, prod(dims) * nsub)
        end
    end
end
