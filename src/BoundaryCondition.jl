"""
    AbstractBoundaryCondition

Abstract supertype for lattice boundary conditions.

This file ships only the minimal vocabulary (`PBC`, `OBC`) required to
exercise the reference lattices (`LineLattice`, `SimpleSquareLattice`)
and the Monte Carlo smoke tests.

The richer per-axis design from `dev/note/04_architecture/03_boundary_and_coordinates`
—`LatticeBoundary{N, A, M}`, `AbstractAxisBC` (`PeriodicAxis` /
`OpenAxis` / `TwistedAxis`), `AbstractBoundaryModifier`—will land in a
subsequent PR. Concrete lattice types here treat `PBC` / `OBC` as
uniform markers that apply to every axis at once.
"""
abstract type AbstractBoundaryCondition end

"""
    PBC()

Periodic boundary condition marker. Applied uniformly to all axes by
the reference lattices.
"""
struct PBC <: AbstractBoundaryCondition end

"""
    OBC()

Open boundary condition marker. Applied uniformly to all axes by the
reference lattices.
"""
struct OBC <: AbstractBoundaryCondition end
