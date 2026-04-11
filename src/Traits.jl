# Lattice trait types.
#
# Only type declarations live here; defaults attached to
# `AbstractLattice` are in `AbstractLattice.jl`.

# ---- Topology ---------------------------------------------------------

"""
    TopologyTrait{Name}

Singleton trait carrying the topology name as a type parameter
(e.g. `TopologyTrait{:Square}`, `TopologyTrait{:Honeycomb}`,
`TopologyTrait{:Penrose}`). Used for dispatch in topology-specific
helpers such as high-symmetry-point tables.
"""
struct TopologyTrait{Name} end

# ---- Periodicity ------------------------------------------------------

"""Periodicity trait marker for Bravais-lattice-like structures."""
struct Periodic end

"""Periodicity trait marker for aperiodic structures (e.g. quasicrystals)."""
struct Aperiodic end

# ---- Reciprocal space support ----------------------------------------

"""
    AbstractReciprocalSupport

Trait indicating what kind of k-space representation a lattice admits.

Subtypes:
- [`HasReciprocal`](@ref): standard Bravais reciprocal lattice.
- [`HasFourierModule`](@ref): dense Fourier module arising from a
  cut-and-project quasicrystal.
- [`NoReciprocal`](@ref): neither of the above (e.g. generic graph).
"""
abstract type AbstractReciprocalSupport end

"""Trait: the lattice has a standard Bravais reciprocal lattice."""
struct HasReciprocal <: AbstractReciprocalSupport end

"""Trait: the lattice has a Fourier module (quasicrystal)."""
struct HasFourierModule <: AbstractReciprocalSupport end

"""Trait: the lattice has no k-space representation."""
struct NoReciprocal <: AbstractReciprocalSupport end

# ---- Size trait -------------------------------------------------------

"""
    AbstractSizeTrait

Trait describing whether a lattice has a finite, infinite, or
finitely-materializable-but-conceptually-infinite extent.

Subtypes:
- [`FiniteSize`](@ref): ordinary finite lattice.
- [`InfiniteSize`](@ref): true infinite lattice (used for
  analytic/spectral work; MC is not applicable).
- [`QuasiInfiniteSize`](@ref): conceptually infinite but materialized
  up to a cutoff (e.g. Penrose radius, Fibonacci depth).
"""
abstract type AbstractSizeTrait end

"""
    FiniteSize{D}(dims::NTuple{D, Int})

Size trait for an ordinary finite lattice with the given per-axis cell
counts.
"""
struct FiniteSize{D} <: AbstractSizeTrait
    dims::NTuple{D, Int}
end

"""
    InfiniteSize()

Size trait for a true infinite lattice. MC cannot run on such a
lattice; spectral / analytic calculations only.
"""
struct InfiniteSize <: AbstractSizeTrait end

"""
    QuasiInfiniteSize{T}(cutoff::T)

Size trait for a conceptually infinite lattice that has been (or will
be) materialized up to a finite cutoff. The cutoff may represent a
radius, a substitution depth, or any other scale parameter native to
the lattice family.
"""
struct QuasiInfiniteSize{T} <: AbstractSizeTrait
    cutoff::T
end
