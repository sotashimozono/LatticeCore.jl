module LatticeCore

using LinearAlgebra
using StaticArrays

# We extend `Base.position` for `AbstractLattice` so downstream code can
# write `position(lat, i)` without name shadowing. `Base.position(io)`
# for IO streams has an unrelated signature.
import Base: position

include("Traits.jl")
include("LatticeElement.jl")
include("SiteType.jl")
include("SiteLayout.jl")
include("AbstractLattice.jl")
include("Bond.jl")
include("BoundaryCondition.jl")
include("Coordinate.jl")
include("Indexing.jl")
include("MomentumLattice.jl")
include("StructureFactor.jl")
include("LazyInfinite.jl")
include("Plot.jl")
include("reference/LineLattice.jl")
include("reference/SimpleSquareLattice.jl")

# ---- AbstractLattice ----
export AbstractLattice
export dimension, num_sites, position, positions, neighbors, boundary
export site_layout, site_type, num_sublattices, sublattice

# ---- Lattice elements ----
export AbstractLatticeElement, VertexCenter, BondCenter, PlaquetteCenter, CellCenter
export element_type

# ---- Site types ----
export AbstractSiteType
export IsingSite, PottsSite, XYSite, HeisenbergSite, EmptySite
export state_type, random_state, zero_state, domain

# ---- Site layouts ----
export AbstractSiteLayout, UniformLayout, SublatticeLayout, ExplicitLayout

# ---- Bond ----
export Bond, bond_center, bonds, neighbor_bonds

# ---- Traits ----
export TopologyTrait, topology
export Periodic, Aperiodic, periodicity
export is_bipartite
export AbstractReciprocalSupport, HasReciprocal, HasFourierModule, NoReciprocal
export reciprocal_support
export AbstractSizeTrait, FiniteSize, InfiniteSize, QuasiInfiniteSize
export size_trait, is_finite

# ---- Boundary conditions ----
export AbstractBoundaryCondition
export AbstractAxisBC, PeriodicAxis, OpenAxis, TwistedAxis
export AbstractBoundaryModifier, NoModifier, SSD
export LatticeBoundary, apply_axis_bc, axis_phase, bond_weight

# ---- Coordinate systems ----
export AbstractCoordinate, RealSpace, LatticeCoord, HigherDimCoord
export to_real, to_lattice, to_hyper

# ---- Indexing ----
export AbstractIndexing, RowMajor, ColMajor, Snake
export site_index, lattice_coord

# ---- Momentum space ----
export AbstractMomentumLattice, PeriodicMomentumLattice
export num_k_points, k_point, reciprocal_basis
export reciprocal_lattice, fourier_module, momentum_lattice
export monkhorst_pack, gamma_centered
export structure_factor
export AcceptanceWindow, HyperReciprocalLattice, BraggPeakSet

# ---- Lazy / infinite ----
export materialize, require_finite

# ---- Plotting (methods live in LatticeCorePlotsExt) ----
export plot_lattice, plot_bonds!, plot_sites!

# ---- Reference lattices ----
export LineLattice, SimpleSquareLattice, basis_vectors

end # module LatticeCore
