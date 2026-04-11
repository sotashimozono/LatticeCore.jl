# Momentum space

This guide is a work in progress. See the concept column
[Reciprocal lattice and Brillouin zone](../concepts/reciprocal_and_brillouin.md)
for the mathematical background.

## AbstractMomentumLattice

k-space is modelled as a subtype of
[`AbstractMomentumLattice{D, T}`](@ref), which itself inherits from
`AbstractLattice{D, T}`. This means the `num_sites`, `position`,
trait, and test-suite infrastructure written for real-space
lattices **also** applies to k-space lattices — k-points are
"sites" and k-vectors are "positions".

## Two concrete subtypes

- [`PeriodicMomentumLattice`](@ref) — eager k-point list built
  from a Bravais reciprocal basis. Constructed via
  [`monkhorst_pack`](@ref) or [`gamma_centered`](@ref).
- [`BraggPeakSet`](@ref) — the quasicrystal counterpart. Same
  `AbstractMomentumLattice` interface, so the same
  [`structure_factor`](@ref) code paths work for both.

## Entry points

- [`reciprocal_lattice(lat)`](@ref) — periodic-only.
- [`fourier_module(lat)`](@ref) — quasicrystal-only (lives in
  `QuasiCrystal.jl`).
- [`momentum_lattice(lat)`](@ref) — trait-dispatched: picks the
  right one via [`reciprocal_support`](@ref).

## Structure factor

[`structure_factor`](@ref) evaluates

```math
S(\mathbf{k}) = \frac{1}{N}\,
    \left| \sum_i s_i \, e^{-i \mathbf{k} \cdot \mathbf{r}_i} \right|^2
```

naively in O(N) per k. FFT-based specialisations can be added
without changing the API.

## TODO

- End-to-end example: compute `S(π, π)` on a 16×16 PBC square for
  the Néel configuration.
- How to sample a k-path (high-symmetry lines) once that lands.
- When to use `monkhorst_pack` vs `gamma_centered`.

## See also

- [Concepts: Reciprocal lattice](../concepts/reciprocal_and_brillouin.md)
- [Reference: Momentum space](../reference/momentum.md)
