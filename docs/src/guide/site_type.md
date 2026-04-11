# Site types and layouts

This guide is a work in progress. See the concept column
[Site types and spin models](../concepts/site_types_and_spins.md)
for the physics.

## Site type ≠ sublattice

LatticeCore separates two ideas that are traditionally conflated:

- **Geometric sublattice** — where inside a unit cell the site
  sits (A, B, C, ...). Stored in `LatticeCoord.sublattice`.
- **Site type** — what physical degree of freedom lives on the
  site (Ising spin, XY rotor, Heisenberg vector, ...). Described
  by [`AbstractSiteType`](@ref) and stored via an
  [`AbstractSiteLayout`](@ref).

These axes are orthogonal: the same geometric sublattice can host
any site type, and the same site type can sit on any sublattice.

## Three layouts

| Layout                         | When to use                                |
| ------------------------------ | ------------------------------------------ |
| [`UniformLayout`](@ref)        | Every site has the same type (most models) |
| [`SublatticeLayout`](@ref)     | Mixed-spin: A = Ising, B = XY, ...         |
| [`ExplicitLayout`](@ref)       | Per-site heterogeneity (quenched disorder) |

`UniformLayout` has zero per-site memory overhead and lets the MC
hot loop constant-fold the site type.

## Built-in site types

[`IsingSite`](@ref), [`PottsSite{Q}`](@ref), [`XYSite`](@ref),
[`HeisenbergSite`](@ref), [`EmptySite`](@ref) — see the
[site type reference](../reference/site_type.md) for
`state_type` / `random_state` / `domain` of each.

## Element centering

The [`element_type`](@ref) trait lets a site type declare that its
DOF lives on a [`BondCenter`](@ref), [`PlaquetteCenter`](@ref), or
[`CellCenter`](@ref) rather than the default
[`VertexCenter`](@ref). This reserves headroom for dimer variables,
gauge links, and flux fields without a breaking rewrite.

## TODO

- Mixed-spin example end-to-end: build a honeycomb with Ising A /
  XY B, then compute `site_type(lat, i)` for each i.
- Custom site type tutorial: add a bond-centered dimer variable.
- Integration with Monte Carlo models in `Lattice2DMonteCarlo.jl`.

## See also

- [Concepts: Site types and spin models](../concepts/site_types_and_spins.md)
- [Reference: Site type](../reference/site_type.md)
- [Reference: Site layout](../reference/site_layout.md)
- [Reference: Lattice element](../reference/lattice_element.md)
