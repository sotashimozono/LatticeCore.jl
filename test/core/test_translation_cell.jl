using LatticeCore
using LatticeCore: CellSite, CellBond
using StaticArrays
using Test

@testset "CellSite / CellBond construction" begin
    s = CellSite((3, -2))
    @test s isa CellSite{2}
    @test s.cell == SVector(3, -2)
    @test s.basis == 1
    @test CellSite((1, 1), 2).basis == 2
    @test CellSite(SVector(0, 0), 1) == CellSite((0, 0), 1)

    b = CellBond(1, 1, (1, 0))
    @test b isa CellBond{2}
    @test b.src == 1 && b.dst == 1
    @test b.offset == SVector(1, 0)
    @test b.type == :nearest
end

@testset "InfiniteSquareLattice — motif / orbits" begin
    inf = InfiniteSquareLattice()

    @test size_trait(inf) == InfiniteSize()
    @test !is_finite(inf)
    @test periodicity(inf) == Periodic()
    @test is_bipartite(inf)
    @test dimension(inf) == 2

    # No finite site count; linear-index API is undefined.
    @test_throws DomainError num_sites(inf)
    @test_throws DomainError position(inf, 1)
    @test_throws DomainError neighbors(inf, 1)

    # Fundamental domain: one site orbit, two bond orbits.
    @test collect(site_orbits(inf)) == [1]
    @test num_basis_sites(inf) == 1
    bo = collect(bond_orbits(inf))
    @test length(bo) == 2
    @test Set(b.offset for b in bo) == Set([SVector(1, 0), SVector(0, 1)])
    @test translation_vectors(inf) == SMatrix{2,2,Float64}(1, 0, 0, 1)
end

@testset "InfiniteSquareLattice — lazy on-demand access" begin
    inf = InfiniteSquareLattice()

    # Positions computed for arbitrary cells without materialising.
    @test cell_position(inf, (0, 0)) == SVector(0.0, 0.0)
    @test cell_position(inf, (3, -2)) == SVector(3.0, -2.0)
    @test cell_position(inf, CellSite((10, 7))) == SVector(10.0, 7.0)

    s = CellSite((0, 0))
    nbrs = neighbors_at(inf, s)
    @test length(nbrs) == 4
    offsets = Set(n.cell for n in nbrs)
    @test offsets == Set(SVector.([(1, 0), (-1, 0), (0, 1), (0, -1)]))
    @test all(n.basis == 1 for n in nbrs)

    # Coordination shell is translation invariant.
    for cell in [(5, 5), (-3, 8), (100, -100)]
        s2 = CellSite(cell)
        shell = Set(n.cell - s2.cell for n in neighbors_at(inf, s2))
        @test shell == Set(SVector.([(1, 0), (-1, 0), (0, 1), (0, -1)]))
    end

    # Incident bonds are anchored at the queried site and self-describing.
    ib = incident_cell_bonds(inf, s)
    @test length(ib) == 4
    @test all(b.src == 1 for b in ib)
    @test Set(b.offset for b in ib) == Set(SVector.([(1, 0), (-1, 0), (0, 1), (0, -1)]))
    # Each incident bond points to the matching neighbour.
    for b in ib
        @test CellSite(s.cell + b.offset, b.dst) in nbrs
    end
end

@testset "materialize bridge (infinite → finite PBC)" begin
    inf = InfiniteSquareLattice()
    fin = materialize(inf; dims=(4, 4))
    @test fin isa SimpleSquareLattice
    @test num_sites(fin) == 16
    @test is_finite(fin)
    @test periodicity(fin) == Periodic()
    # A materialised interior site has the same coordination as the motif.
    @test length(neighbors(fin, 6)) == 4
end

@testset "uniform orbit UI across finite and infinite" begin
    fin = SimpleSquareLattice(4, 4, PeriodicAxis())
    inf = InfiniteSquareLattice()

    # Same fundamental domain regardless of finiteness / BC.
    @test collect(site_orbits(fin)) == collect(site_orbits(inf))
    @test num_basis_sites(fin) == num_basis_sites(inf)

    fin_bonds = Set((b.src, b.dst, b.offset) for b in bond_orbits(fin))
    inf_bonds = Set((b.src, b.dst, b.offset) for b in bond_orbits(inf))
    @test fin_bonds == inf_bonds
    @test translation_vectors(fin) == translation_vectors(inf)

    # The motif reproduces the finite lattice's own bond displacements:
    # a PBC interior site steps to the four unit neighbours.
    inf_shell = Set(b.offset for b in incident_cell_bonds(inf, CellSite((0, 0))))
    @test inf_shell == Set(SVector.([(1, 0), (-1, 0), (0, 1), (0, -1)]))
end
