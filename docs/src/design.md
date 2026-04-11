# Design notes

This page summarises the architectural decisions that shape
LatticeCore. The long-form notes — chapter-by-chapter rationale,
alternatives considered, cross-references between chapters — live
outside the package in the parent workspace's
`dev/note/04_architecture/` tree.

## Five orthogonal axes

LatticeCore models every lattice along five independent axes. The
cross-references point to the chapter of `dev/note/04_architecture/`
where each axis is designed.

| Axis                    | Lives in                                                    | Chapter |
| ----------------------- | ----------------------------------------------------------- | ------- |
| Topology / graph        | `Traits.jl` + concrete lattice's `neighbors`                | 02      |
| Boundary condition      | [`LatticeBoundary`](@ref), [`AbstractAxisBC`](@ref)         | 03      |
| Coordinate system       | [`AbstractCoordinate`](@ref) + [`AbstractIndexing`](@ref)   | 03      |
| Geometric sublattice    | `LatticeCoord.sublattice` + `sublattice(lat, i)`            | 04      |
| Physical site type      | [`AbstractSiteType`](@ref) + [`AbstractSiteLayout`](@ref)   | 04      |

Momentum space (chapter 05) and the size trait hierarchy /
materialisation pattern (chapter 06) sit on top of these five axes
and reuse the same interface machinery.

## Key rules

1. **Two type parameters only.** `AbstractLattice{D, T}`. The rest
   is field storage and dispatch, so the type system stays small.
2. **No `isa` branching in hot paths.** Every variant point goes
   through multiple dispatch or a Holy trait. This keeps extensions
   additive: a new axis BC is a new subtype plus a new method, not
   an edit to an `if` chain.
3. **Sublattice and site type are different concepts.** The
   geometric sublattice is where a site sits inside the unit cell;
   the site type is what lives on it. See
   [Lattices and unit cells](concepts/lattice_and_unit_cells.md).
4. **Per-axis boundary conditions.** Cylinders (PBC × OBC) and
   twisted periodic lattices drop out of the same type system as
   uniform PBC.
5. **Infinite structures are expressed abstractly and materialised
   at a cutoff.** The same pattern applies in real space
   ([`materialize`](@ref)) and in reciprocal space
   ([`BraggPeakSet`](@ref) as a finite snapshot of
   [`HyperReciprocalLattice`](@ref)).
6. **Lazy-friendly interface.** `position(lat, i)` and
   `neighbors(lat, i)` are pinpoint; `positions(lat)` and
   `bonds(lat)` are iterators. No method requires enumerating
   every site up front.
7. **Breaking changes are cheap in 0.x.** The whole stack is
   pre-1.0; the `VersionCheck` workflow enforces a version bump on
   every PR, and tags track every merged change.

## Where the notes live

The authoritative design chapters are under `dev/note/04_architecture/`
in the parent workspace:

- `01_package_layout/` — why LatticeCore is a separate repo
- `02_lattice_interface/` — the `AbstractLattice` contract
- `03_boundary_and_coordinates/` — per-axis BC and coordinate system
- `04_site_type/` — site types, layouts, and element centering
- `05_momentum_space/` — reciprocal lattices, Bragg peaks, structure
  factor
- `06_lazy_infinite/` — size traits and materialisation
- `07_mc_layer/` — Monte Carlo consumer layer (lives in
  `Lattice2DMonteCarlo.jl`)

They are checked into the parent workspace rather than into
LatticeCore itself so that cross-package architectural changes can
be discussed in one place.
