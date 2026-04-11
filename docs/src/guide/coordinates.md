# Coordinate systems and indexing

This guide is a work in progress. See the
[coordinate API reference](../reference/coordinates.md) and the
[indexing reference](../reference/indexing.md) for the exact
signatures.

## Three coordinate kinds

LatticeCore distinguishes three ways of naming a lattice point:

- [`RealSpace{D, T}`](@ref) — Cartesian coordinates, stored as an
  `SVector{D, T}`.
- [`LatticeCoord{D}`](@ref) — an integer cell index plus a 1-based
  geometric sublattice id.
- [`HigherDimCoord{DPhys, DHyper, T}`](@ref) — a point in the host
  lattice of a cut-and-project quasicrystal.

They talk to each other through [`to_real`](@ref),
[`to_lattice`](@ref), and [`to_hyper`](@ref), all of which dispatch
on the concrete lattice so that the conversion uses the lattice's
basis and (if applicable) projection.

## Indexing is orthogonal to coordinates

Once you have a [`LatticeCoord`](@ref), you still need a rule for
turning it into a 1-based site index — that is the job of the
[`AbstractIndexing`](@ref) hierarchy. LatticeCore ships
[`RowMajor`](@ref), [`ColMajor`](@ref), and [`Snake`](@ref) with
exhaustive round-trip tests.

Indexing is deliberately **decoupled** from coordinates: a lattice
can switch its indexing without changing how it talks about
positions.

## TODO

- Worked example: writing a custom indexing strategy (e.g. Hilbert
  curve) and verifying the round-trip with `site_index` /
  `lattice_coord`.
- When the sublattice-innermost convention is the wrong choice and
  how to override it.
- HigherDimCoord and the cut-and-project projection matrices for
  Penrose / Ammann–Beenker / Fibonacci.

## See also

- [Concepts: Lattices and unit cells](../concepts/lattice_and_unit_cells.md)
- [Reference: Coordinates](../reference/coordinates.md)
- [Reference: Indexing](../reference/indexing.md)
