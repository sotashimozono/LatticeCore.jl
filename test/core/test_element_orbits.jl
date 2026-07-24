using LatticeCore
using LatticeCore: CellBond
using StaticArrays
using Test

@testset "element_orbits routes centring to the motif orbits" begin
    inf = InfiniteSquareLattice()

    @test collect(element_orbits(inf, VertexCenter())) == collect(site_orbits(inf))
    @test collect(element_orbits(inf, BondCenter())) == collect(bond_orbits(inf))

    pl = collect(element_orbits(inf, PlaquetteCenter()))
    @test length(pl) == 1
    @test pl[1].type == :square
    @test length(pl[1].corners) == 4
end

@testset "element centring is accessible on the infinite lattice (no throw)" begin
    inf = InfiniteSquareLattice()
    # The enumeration API throws (no linear index); the orbit API must not.
    @test_throws DomainError elements(inf, VertexCenter())         # 1:num_sites
    @test element_orbits(inf, VertexCenter()) == 1:1
    @test length(collect(element_orbits(inf, BondCenter()))) == 2
    @test length(collect(element_orbits(inf, PlaquetteCenter()))) == 1
end

@testset "orbit representative positions are geometrically consistent" begin
    inf = InfiniteSquareLattice()

    # Vertex: the single basis site sits at the cell origin.
    @test element_orbit_position(inf, VertexCenter(), 1) == SVector(0.0, 0.0)

    # Bond: midpoint of the two endpoints.
    bx, by = collect(bond_orbits(inf))
    @test element_orbit_position(inf, BondCenter(), bx) == SVector(0.5, 0.0)
    @test element_orbit_position(inf, BondCenter(), by) == SVector(0.0, 0.5)
    # Independent check: the midpoint equals mean of the two endpoint
    # cell positions, reconstructed directly.
    for cb in (bx, by)
        p =
            (
                cell_position(inf, (0, 0), cb.src) +
                cell_position(inf, Tuple(cb.offset), cb.dst)
            ) / 2
        @test element_orbit_position(inf, BondCenter(), cb) == p
    end

    # Plaquette: centroid of the unit square is its centre.
    sq = only(element_orbits(inf, PlaquetteCenter()))
    @test element_orbit_position(inf, PlaquetteCenter(), sq) == SVector(0.5, 0.5)
    # Independent check: centroid == mean of the four corner positions.
    corners = [cell_position(inf, (dx, dy), sub) for (sub, dx, dy) in sq.corners]
    @test element_orbit_position(inf, PlaquetteCenter(), sq) == sum(corners) / 4
end

@testset "finite and infinite agree on every centring" begin
    fin = SimpleSquareLattice(4, 4, PeriodicAxis())
    inf = InfiniteSquareLattice()

    @test collect(element_orbits(fin, VertexCenter())) ==
        collect(element_orbits(inf, VertexCenter()))
    fkey(cb) = (cb.src, cb.dst, cb.offset, cb.type)
    @test Set(fkey.(element_orbits(fin, BondCenter()))) ==
        Set(fkey.(element_orbits(inf, BondCenter())))

    fpl = only(element_orbits(fin, PlaquetteCenter()))
    ipl = only(element_orbits(inf, PlaquetteCenter()))
    @test fpl.corners == ipl.corners && fpl.type == ipl.type

    @test element_orbit_position(fin, VertexCenter(), 1) ==
        element_orbit_position(inf, VertexCenter(), 1)
    for cb in element_orbits(inf, BondCenter())
        @test element_orbit_position(fin, BondCenter(), cb) ==
            element_orbit_position(inf, BondCenter(), cb)
    end
    @test element_orbit_position(fin, PlaquetteCenter(), fpl) ==
        element_orbit_position(inf, PlaquetteCenter(), ipl)
end

@testset "plaquette_orbits defaults to empty for lattices without faces" begin
    line = LineLattice(5)
    @test plaquette_orbits(line) == ()
    @test element_orbits(line, PlaquetteCenter()) == ()
end
