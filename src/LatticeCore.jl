module LatticeCore

using LinearAlgebra
using StaticArrays

# We extend `Base.position` for `AbstractLattice` so downstream code can
# write `position(lat, i)` without name shadowing. `Base.position(io)`
# for IO streams has an unrelated signature.
import Base: position

include("Traits.jl")
include("AbstractLattice.jl")
include("Bond.jl")
include("BoundaryCondition.jl")
include("reference/LineLattice.jl")
include("reference/SimpleSquareLattice.jl")

# ---- AbstractLattice ----
export AbstractLattice
export dimension, num_sites, position, positions, neighbors, boundary

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
export AbstractBoundaryCondition, PBC, OBC

# ---- Reference lattices ----
export LineLattice, SimpleSquareLattice

end # module LatticeCore
