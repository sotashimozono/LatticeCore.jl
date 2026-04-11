"""
    AbstractLattice{D, T}

Abstract supertype for lattice types. Type parameters:
- `D`: spatial dimension (integer literal)
- `T`: numeric type for real-space positions (typically `Float64`)

# Required interface (concrete subtypes must implement)
- `position(lat, i)::SVector{D, T}`
- `neighbors(lat, i)` — iterable of neighbor site indices
- `boundary(lat)` — returns an `AbstractBoundaryCondition`
- `size_trait(lat)::AbstractSizeTrait`
- `num_sites(lat)::Int` — required for `FiniteSize` / `QuasiInfiniteSize`
  lattices; `InfiniteSize` lattices should throw `DomainError`.

# Optional interface (defaults provided)
- `positions(lat)` — defaults to an iterator built from `position`
- `bonds(lat)` — defaults to an iterator built from `neighbors`
- `neighbor_bonds(lat, i)` — defaults to an iterator built from `neighbors`
- `topology(lat)` — `TopologyTrait{:unknown}()`
- `periodicity(lat)` — `Aperiodic()`
- `is_bipartite(lat)` — `false`
- `reciprocal_support(lat)` — `NoReciprocal()`
- `is_finite(lat)` — derived from `size_trait`

See `dev/note/04_architecture/02_lattice_interface/README.md` for the
design rationale.
"""
abstract type AbstractLattice{D,T} end

# ---- Dimension helpers (fully determined by type parameters) ---------

"""
    dimension(lat) → Int

Spatial dimension `D` of the lattice.
"""
dimension(::AbstractLattice{D,T}) where {D,T} = D

Base.eltype(::Type{<:AbstractLattice{D,T}}) where {D,T} = T
Base.eltype(lat::AbstractLattice) = eltype(typeof(lat))

# ---- Required interface (generic functions, methods added per type) --

"""
    num_sites(lat::AbstractLattice)::Int

Number of sites in the lattice. Must be implemented by finite concrete
types. `InfiniteSize` lattices should throw `DomainError`.
"""
function num_sites end

"""
    position(lat::AbstractLattice{D, T}, i::Int)::SVector{D, T}

Real-space position of site `i`. This extends `Base.position` so that
downstream code can write `position(lat, i)` directly after
`using LatticeCore`.
"""
position(lat::AbstractLattice, i::Int) = throw(MethodError(position, (lat, i)))

"""
    neighbors(lat::AbstractLattice, i::Int)

Indices of sites adjacent to site `i` (nearest neighbors by default).
Concrete types may additionally define `neighbors(lat, i; shell)` or
`neighbors(lat, i, shell::Int)` to support higher shells.
"""
function neighbors end

"""
    boundary(lat::AbstractLattice)

The lattice's boundary condition (subtype of
`AbstractBoundaryCondition`, defined in `BoundaryCondition.jl`).
"""
function boundary end

"""
    size_trait(lat::AbstractLattice) → AbstractSizeTrait

Size trait describing the lattice's extent. Must be implemented by
concrete types.
"""
function size_trait end

# ---- Optional interface with defaults --------------------------------

"""
    positions(lat::AbstractLattice)

Iterator over all positions. Default implementation lazily constructs
from `num_sites` and `position`; works only for finite lattices.
Concrete types may override for efficiency.
"""
function positions(lat::AbstractLattice)
    return (position(lat, i) for i in 1:num_sites(lat))
end

# ---- Trait defaults --------------------------------------------------

"""Default topology trait: `TopologyTrait{:unknown}()`."""
topology(::AbstractLattice) = TopologyTrait{:unknown}()

"""Default periodicity: `Aperiodic()`."""
periodicity(::AbstractLattice) = Aperiodic()

"""Default bipartite flag: `false`."""
is_bipartite(::AbstractLattice) = false

"""Default reciprocal support: `NoReciprocal()`."""
reciprocal_support(::AbstractLattice) = NoReciprocal()

"""
    is_finite(lat::AbstractLattice) → Bool

`true` if `size_trait(lat)` is a [`FiniteSize`](@ref). Monte Carlo
algorithms should guard against non-finite lattices with this
predicate.
"""
is_finite(lat::AbstractLattice) = size_trait(lat) isa FiniteSize
