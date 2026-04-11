# Site types and spin models

This column is still being written. The short version:

## What is a site type?

A **site type** tells LatticeCore what mathematical object lives on
a single lattice site and how to sample it uniformly at random.
Concretely, [`AbstractSiteType`](@ref) requires

- [`state_type`](@ref) — the Julia type that stores one state
- [`random_state`](@ref) — a uniform sampler

and optionally

- [`zero_state`](@ref), [`domain`](@ref), and
  [`element_type`](@ref) for lifting the DOF onto bonds / plaquettes
  / cells.

## The classical spin ladder

The archetypal site types form a ladder of increasing continuous
symmetry:

| Site type                | State space                | Symmetry group | Models            |
| ------------------------ | -------------------------- | -------------- | ----------------- |
| [`IsingSite`](@ref)      | ``\{-1, +1\}``             | ``\mathbb{Z}_2`` | Ising, RBIM     |
| [`PottsSite{Q}`](@ref)   | ``\{1, \dots, Q\}``        | ``S_Q``         | Potts, clock     |
| [`XYSite`](@ref)         | circle ``S^1``             | ``U(1) = SO(2)``  | XY, Villain   |
| [`HeisenbergSite`](@ref) | sphere ``S^2 \subset \mathbb{R}^3`` | ``SO(3)`` | Heisenberg |
| [`EmptySite`](@ref)      | singleton                  | trivial         | vacancies       |

The pattern is standard condensed-matter: as you move down the
table, the Mermin–Wagner theorem kicks in (no long-range order in
2D for ``d \geq 2`` continuous groups at finite ``T``), the energy
landscape becomes smoother, and the numerically appropriate Monte
Carlo update changes (single-flip → cluster → overrelaxation,
cluster → worm).

## Why this lives in LatticeCore

A Monte Carlo runtime like `Lattice2DMonteCarlo.jl` needs to
dispatch its local Hamiltonian on the site types at either end of a
bond. Making `AbstractSiteType` part of LatticeCore — rather than
the MC layer — means that *every* concrete lattice in the stack can
describe what it carries, regardless of which MC runtime later
consumes it.

## TODO

- Worked example: how the Hamiltonian dispatch actually works when
  A and B sublattices carry different site types (`SublatticeLayout`).
- Heisenberg uniform sphere sampler (Marsaglia / Archimedes) and
  why `HeisenbergSite` uses `SVector{3, Float64}`.
- Mermin–Wagner in one paragraph and how it motivates the
  Ising → Potts → XY → Heisenberg escalation.
- How to build a custom site type for a bond-centered DOF
  (dimer model, gauge link) via [`element_type`](@ref).
