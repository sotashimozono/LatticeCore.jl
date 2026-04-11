# Boundary conditions

This column is still being written. In the meantime, the physics
that drives LatticeCore's per-axis [`LatticeBoundary`](@ref)
abstraction:

- **Periodic (PBC)**: the lattice is wrapped onto a torus. It is
  the gentlest finite-size approximation — every site looks like a
  bulk site — and it lets the reciprocal lattice remain discrete.
- **Open (OBC)**: the sample has real surfaces. Edge / corner sites
  have fewer neighbours, and the discrete reciprocal lattice no
  longer strictly applies.
- **Cylinder / strip**: a mix of PBC and OBC along different axes.
  Ubiquitous in DMRG and tensor-network simulations because it
  tames edge effects while keeping a manageable transverse size.
- **Twisted PBC**: periodic with a phase factor
  ``e^{i\theta}`` attached to bonds that cross the boundary. This
  is how a lattice model sees a magnetic flux, a twist angle, or a
  supercurrent. The topology is unchanged; only the phase
  accumulates.
- **Sine-square deformation (SSD)**: a non-topological *modifier*
  that multiplies bond energies by a smooth envelope
  ``\sin^2(\pi x / L)``. SSD is topologically open but reproduces
  periodic-like ground states surprisingly well, which is useful
  when you want PBC physics without the PBC integration pain.

## Why per-axis?

Classical physics codes usually hardcode "PBC" or "OBC" as a global
flag. That is wrong as soon as you want a cylinder — the axis
choices stop being symmetric. LatticeCore follows the convention
used by DMRG / tensor-network libraries and stores **one
[`AbstractAxisBC`](@ref) per axis**, plus one
[`AbstractBoundaryModifier`](@ref), in a composite
[`LatticeBoundary`](@ref). Adding a new BC kind means adding a new
`AbstractAxisBC` subtype and a method for
[`apply_axis_bc`](@ref); no core code changes.

## TODO

- Diagram: PBC torus, OBC square, cylinder, Möbius.
- Worked example: dispersion of a 1D tight-binding chain under
  `PeriodicAxis`, `OpenAxis`, and `TwistedAxis`.
- Finite-size scaling and the role of BC in extracting
  bulk-limit quantities.
- Where SSD fits in the modifier slot of [`LatticeBoundary`](@ref).
