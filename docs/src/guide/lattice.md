# Lattice interface

This guide walks through the `AbstractLattice` interface — the
minimum a concrete lattice must implement, the trait-based extension
points, and the default methods that come for free.

For the physics and mathematical background, see the concept column
[Lattices and unit cells](../concepts/lattice_and_unit_cells.md).

## The abstract type

```julia
abstract type AbstractLattice{D, T} end
```

The two type parameters are:

- `D::Int` — the physical dimension (1, 2, 3, ...)
- `T<:Real` — the numeric type used to store positions (typically
  `Float64`)

**Nothing else** lives in the type parameters. Boundary conditions,
topology, indexing, and site types are all stored as fields on the
concrete subtype, and dispatched through multiple dispatch and
trait objects.

## Required interface

A concrete subtype must implement these:

| Method                 | Returns                                |
| ---------------------- | -------------------------------------- |
| [`num_sites`](@ref)    | `Int`                                  |
| [`position`](@ref)     | `SVector{D, T}`                        |
| [`neighbors`](@ref)    | `AbstractVector{Int}` (pinpoint)       |
| [`boundary`](@ref)     | [`LatticeBoundary`](@ref)              |
| [`site_layout`](@ref)  | [`AbstractSiteLayout`](@ref)           |
| [`size_trait`](@ref)   | [`AbstractSizeTrait`](@ref)            |

`num_sites` may be undefined — or may throw — when the lattice's
`size_trait` is [`InfiniteSize`](@ref).

## Free defaults

Once the required methods are in place, LatticeCore derives:

- [`dimension`](@ref) and `Base.eltype` from the type parameters
- [`positions`](@ref) and [`bonds`](@ref) as lazy iterators from
  `position` / `neighbors`
- [`neighbor_bonds`](@ref) from `neighbors`
- [`is_finite`](@ref), [`is_bipartite`](@ref),
  [`reciprocal_support`](@ref), [`periodicity`](@ref),
  [`topology`](@ref) as safe defaults that can be overridden
- [`site_type(lat, i)`](@ref site_type) by delegating to
  `site_layout`

## Trait-based extension

Extension points use value types rather than flags, so dispatch
stays compile-time. The main traits are:

- [`TopologyTrait`](@ref) — `{:square}`, `{:honeycomb}`, etc.
- [`Periodic`](@ref) / [`Aperiodic`](@ref) — is the underlying
  crystal periodic?
- [`AbstractReciprocalSupport`](@ref) — does this lattice admit a
  Bravais reciprocal ([`HasReciprocal`](@ref)), a dense Fourier
  module ([`HasFourierModule`](@ref)), or neither
  ([`NoReciprocal`](@ref))?
- [`AbstractSizeTrait`](@ref) — [`FiniteSize`](@ref),
  [`InfiniteSize`](@ref), or [`QuasiInfiniteSize`](@ref).

## API reference

- [Lattice](../reference/lattice.md)
- [Bond](../reference/bond.md)
- [Traits](../reference/traits.md)
- [Reference lattices](../reference/reference_lattices.md)
