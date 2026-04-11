# LatticeCore.jl

LatticeCore is the **abstract interface layer** for a family of Julia
packages that simulate physical systems on geometric lattices. It
defines the types and trait vocabulary that every lattice — periodic
or aperiodic, finite or conceptually infinite — must implement, and
then ships a pair of minimal reference implementations
([`LineLattice`](@ref), [`SimpleSquareLattice`](@ref)) plus a
momentum-space layer so the interface is verifiable end-to-end
without any other dependency.

## What LatticeCore is

- A single abstract type hierarchy rooted at
  `AbstractLattice{D, T}`, parameterised **only** on the physical
  dimension `D` and the numeric type `T`. Everything else — boundary
  conditions, topology, indexing, site types — is a field, composed
  via multiple dispatch and trait objects.
- A **trait-based extension model** (`TopologyTrait`, `Periodic` /
  `Aperiodic`, `HasReciprocal` / `HasFourierModule` / `NoReciprocal`,
  the size trait hierarchy) so new lattice kinds can plug in without
  touching core types.
- A **per-axis boundary condition** abstraction
  ([`LatticeBoundary`](@ref)) so cylinders (PBC × OBC) and
  flux-twisted lattices are as easy to describe as uniform PBC.
- A clean split between **geometric sublattice** (A/B on a honeycomb,
  A/B/C on kagome) and **site type** (the physical DOF living on the
  site: Ising spin, XY rotor, Heisenberg vector, vacancy).
- A **momentum-space layer** ([`AbstractMomentumLattice`](@ref),
  [`PeriodicMomentumLattice`](@ref), [`structure_factor`](@ref))
  plus the type skeletons [`HyperReciprocalLattice`](@ref) /
  [`BraggPeakSet`](@ref) so quasicrystal Fourier analysis can ride
  the same observer codepaths as periodic lattices.
- A **lazy-friendly interface**: `position(lat, i)` and
  `neighbors(lat, i)` are pinpoint, `positions(lat)` and `bonds(lat)`
  are iterators. Infinite structures can be expressed abstractly and
  lowered to finite lattices through [`materialize`](@ref).

## What LatticeCore is **not**

- A production-grade 2D lattice catalogue (triangular, honeycomb,
  kagome, Lieb, Shastry–Sutherland, ...). Those live in
  `Lattice2D.jl` and are built on top of LatticeCore.
- A quasicrystal generator. Cut-and-project Penrose / Ammann–Beenker
  / Fibonacci, their acceptance windows, and the Bragg-peak enumerator
  live in `QuasiCrystal.jl`.
- A Monte Carlo runtime. Classical MC models and algorithms live in
  `Lattice2DMonteCarlo.jl`; LatticeCore only supplies the site-type
  system and the [`require_finite`](@ref) guard they rely on.
- A visualisation layer. Plots / Makie integration lives in package
  extensions on top of the concrete-lattice packages, not here.

## Package layout

```text
LatticeCore.jl            (this package, abstract interface)
      │
      ├── Lattice2D.jl        (periodic 2D catalogue)
      │        │
      ├── QuasiCrystal.jl     (cut-and-project quasicrystals)
      │        │
      └── ────┬─┘
              ▼
     Lattice2DMonteCarlo.jl   (classical MC runtime)
```

## Design principles

1. **One abstract type, minimal parameters.** `AbstractLattice{D, T}`.
   BC / site type / indexing live in fields.
2. **No `isa` branching.** Every variant goes through multiple
   dispatch or a Holy trait.
3. **Separation of concerns.** Topology (which bonds exist), boundary
   (how they wrap), modifier (how they reweight), geometric
   sublattice, and physical site type are five independent axes.
4. **Lazy-friendly I/F.** Iterators everywhere; no method requires a
   full site enumeration.
5. **Breaking changes are cheap** in the 0.x era. The whole stack is
   pre-1.0 and evolves with the architecture notes in the parent
   repo.

## Further reading

- [Getting Started](getting_started.md) — install, construct a
  lattice, measure it.
- [Guide](guide/lattice.md) — concept-by-concept walkthrough.
- [API Reference](reference/lattice.md) — generated from docstrings.
- [Design Notes](design.md) — summary of the key architectural
  decisions and pointers to the long-form notes.
