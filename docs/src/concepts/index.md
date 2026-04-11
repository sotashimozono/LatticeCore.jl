# Concepts

This section is a set of **columns**: short, self-contained reads
that introduce the concepts behind LatticeCore's design. They are
written for someone approaching lattice simulation for the first
time — or from the software-engineering side — and want to understand
*why* the API is shaped the way it is, not just how to call it.

The goal is not to re-derive condensed-matter physics; it is to make
the vocabulary LatticeCore borrows from physics legible, and to show
where each concept surfaces in the code.

## Where to read first

If you are new to lattices and are here to simulate something, read
in order:

1. [Lattices and unit cells](lattice_and_unit_cells.md) — what a
   Bravais lattice is, how unit cells and sublattices work, and why
   LatticeCore treats the *geometric sublattice* and the *physical
   site type* as independent axes.
2. [Reciprocal lattice and Brillouin zone](reciprocal_and_brillouin.md)
   — the dual picture in momentum space. This is the mathematical
   home of `reciprocal_lattice`, Monkhorst–Pack meshes, and the
   structure factor.
3. [Boundary conditions](boundary_conditions.md) — periodic vs open,
   cylinders, and twisted boundaries, and how each changes the
   physics of finite-size simulations.
4. [Site types and spin models](site_types_and_spins.md) — Ising,
   Potts, XY, Heisenberg; what the distinction *means* physically
   and how LatticeCore represents it.
5. [Quasiperiodic order](quasiperiodic_order.md) — cut-and-project
   quasicrystals, Fourier modules, and why the reciprocal lattice
   becomes dense.

## Relation to the rest of the docs

| Concept column            | API guide                             | Reference                             |
| ------------------------- | ------------------------------------- | ------------------------------------- |
| Lattices and unit cells   | [Lattice interface](../guide/lattice.md)     | [`AbstractLattice`](@ref)            |
| Reciprocal / Brillouin    | [Momentum space](../guide/momentum.md)       | [`reciprocal_lattice`](@ref)         |
| Boundary conditions       | [Boundary conditions](../guide/boundary.md)  | [`LatticeBoundary`](@ref)            |
| Site types and spins      | [Site types](../guide/site_type.md)          | [`AbstractSiteType`](@ref)           |
| Quasiperiodic order       | [Lazy / infinite](../guide/lazy_infinite.md) | [`HyperReciprocalLattice`](@ref)     |
