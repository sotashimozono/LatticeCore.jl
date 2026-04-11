# Lazy / infinite lattices

This guide is a work in progress. See the concept column
[Quasiperiodic order](../concepts/quasiperiodic_order.md) for the
physical motivation.

## Size traits

LatticeCore describes the extent of a lattice through an
[`AbstractSizeTrait`](@ref) value attached to every concrete lattice
via [`size_trait`](@ref):

- [`FiniteSize`](@ref) — ordinary finite lattice. MC is safe.
- [`InfiniteSize`](@ref) — conceptually infinite. MC is *not* safe;
  use the lattice only for spectral / analytic work.
- [`QuasiInfiniteSize`](@ref) — conceptually infinite but expected
  to be materialised at a cutoff before use.

[`is_finite(lat)`](@ref) is a one-liner derived from this trait and
is the predicate Monte Carlo algorithms should guard on.

## The infinite-abstract → finite-materialisation pattern

LatticeCore ships two small hooks for this pattern:

- [`materialize(abstract; kwargs...)`](@ref materialize) — a generic
  function with **no typed supertype**. A package that defines an
  "infinite abstract" type (e.g. `InfiniteFibonacci`) adds a
  method for it that returns a `FiniteSize` concrete lattice.
- [`require_finite(lat)`](@ref require_finite) — an assertion guard
  that MC entry points should call:

```julia
function run!(rng, state, lat::AbstractLattice, model, alg; kwargs...)
    require_finite(lat)
    # ... safe from here ...
end
```

This pattern mirrors the momentum-space counterpart:
[`HyperReciprocalLattice`](@ref) is the "infinite abstract" that
becomes a finite [`BraggPeakSet`](@ref) at a cutoff.

## TODO

- Example: a custom `InfiniteChain` type + matching `materialize`
  that produces a `LineLattice` at a chosen depth.
- When real on-demand laziness (as opposed to materialisation)
  makes sense, and the open design questions around it.

## See also

- [Concepts: Quasiperiodic order](../concepts/quasiperiodic_order.md)
- [Reference: Lazy / infinite](../reference/lazy_infinite.md)
