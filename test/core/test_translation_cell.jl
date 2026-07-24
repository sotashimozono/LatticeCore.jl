using LatticeCore
using LatticeCore: CellSite, CellBond
using StaticArrays
using Test

# A minimal two-basis 1D lattice used to exercise the asymmetric-motif
# (`src != dst`) paths of `incident_cell_bonds` that the single-basis
# square lattices cannot distinguish. Basis site A sits at 0.0, B at
# 0.5; A–B are bonded within a cell and B–A across the +1 cell, so the
# chain is …A B A B… with coordination 2.
struct TwoBasisChain <: LatticeCore.AbstractLattice{1,Float64} end
LatticeCore.translation_vectors(::TwoBasisChain) = SMatrix{1,1,Float64}(1.0)
LatticeCore.num_basis_sites(::TwoBasisChain) = 2
function LatticeCore.basis_position(::TwoBasisChain, b::Int)
    b == 1 && return SVector(0.0)
    b == 2 && return SVector(0.5)
    return throw(ArgumentError("bad basis $b"))
end
function LatticeCore.cell_bonds(::TwoBasisChain)
    return (CellBond(1, 2, (0,), :nearest), CellBond(2, 1, (1,), :nearest))
end

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

@testset "basis_position default and error path" begin
    inf = InfiniteSquareLattice()
    @test basis_position(inf, 1) == SVector(0.0, 0.0)
    # The single-basis default rejects any other basis index.
    @test_throws ArgumentError basis_position(inf, 2)
end

@testset "out-of-range basis is rejected, not silently empty" begin
    inf = InfiniteSquareLattice()
    # A malformed CellSite must fail loudly rather than return an empty
    # (wrong) coordination shell.
    @test_throws ArgumentError incident_cell_bonds(inf, CellSite((0, 0), 2))
    @test_throws ArgumentError neighbors_at(inf, CellSite((0, 0), 2))
    @test_throws ArgumentError incident_cell_bonds(inf, CellSite((0, 0), 0))
    # The valid basis still works.
    @test length(neighbors_at(inf, CellSite((0, 0), 1))) == 4
end

@testset "materialize propagates a custom site layout" begin
    layout = UniformLayout(XYSite())
    inf = InfiniteSquareLattice(; layout=layout)
    @test site_layout(inf) == layout
    @test site_type(inf, 1) == XYSite()

    fin = materialize(inf; dims=(3, 3))
    @test site_layout(fin) == layout
    @test site_type(fin, 1) == XYSite()
end

@testset "asymmetric (src != dst) motif — reverse-orientation branch" begin
    chain = TwoBasisChain()
    @test num_basis_sites(chain) == 2
    @test collect(site_orbits(chain)) == [1, 2]

    # Positions within a cell, and across cells, computed lazily.
    @test cell_position(chain, (0,), 1) == SVector(0.0)
    @test cell_position(chain, (0,), 2) == SVector(0.5)
    @test cell_position(chain, (2,), 2) == SVector(2.5)

    # Basis A (=1) at cell 0: forward branch of bond (1,2,(0,)) → B in
    # this cell; reverse branch of bond (2,1,(1,)) → B in the −1 cell.
    ibA = incident_cell_bonds(chain, CellSite((0,), 1))
    @test all(b.src == 1 for b in ibA)
    @test Set((b.dst, b.offset) for b in ibA) == Set([(2, SVector(0)), (2, SVector(-1))])
    nbrA = neighbors_at(chain, CellSite((0,), 1))
    @test Set((n.cell, n.basis) for n in nbrA) == Set([(SVector(0), 2), (SVector(-1), 2)])

    # Basis B (=2) at cell 0: reverse branch of bond (1,2,(0,)) → A in
    # this cell (offset negated to (0,)); forward branch of bond
    # (2,1,(1,)) → A in the +1 cell. This is the path the single-basis
    # square motif cannot exercise.
    ibB = incident_cell_bonds(chain, CellSite((0,), 2))
    @test all(b.src == 2 for b in ibB)
    @test Set((b.dst, b.offset) for b in ibB) == Set([(1, SVector(0)), (1, SVector(1))])
    nbrB = neighbors_at(chain, CellSite((0,), 2))
    @test Set((n.cell, n.basis) for n in nbrB) == Set([(SVector(0), 1), (SVector(1), 1)])

    # Reciprocity: A at cell 0 lists B at cell −1 as a neighbour, and B
    # at cell −1 must list A at cell 0 back.
    nbrB_m1 = neighbors_at(chain, CellSite((-1,), 2))
    @test (SVector(0), 1) in Set((n.cell, n.basis) for n in nbrB_m1)
end
