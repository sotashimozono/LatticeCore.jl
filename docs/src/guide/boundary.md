# Boundary conditions

This guide is a work in progress. For now, see the concept column
[Boundary conditions](../concepts/boundary_conditions.md) for the
physics, and the [boundary API reference](../reference/boundary.md)
for the exact signatures.

## The picture

A boundary condition in LatticeCore is a composite object:

```julia
LatticeBoundary(
    axes     = (PeriodicAxis(), OpenAxis()),   # one AbstractAxisBC per axis
    modifier = NoModifier(),                   # non-topological reweighting
)
```

The **axes** decide which candidate bonds exist — they wrap, clamp,
or attach a phase to a bond that would otherwise leave the sample.
The **modifier** decides how existing bonds are weighted (e.g.
sine-square deformation). Splitting the two concerns this way means
a cylinder and an SSD modification can be composed without any core
code change.

## What LatticeCore ships today

- Axis BCs: [`PeriodicAxis`](@ref), [`OpenAxis`](@ref),
  [`TwistedAxis`](@ref).
- Modifiers: [`NoModifier`](@ref), [`SSD`](@ref) (type only; the
  weighting rule will land with the MC layer).
- Hooks: [`apply_axis_bc`](@ref), [`axis_phase`](@ref),
  [`bond_weight`](@ref).

## TODO

- Worked cylinder example on a 3×4 square lattice, with a picture of
  the neighbour graph.
- How reference lattices use `apply_axis_bc` inside `neighbors`.
- Adding a new axis BC subtype from scratch.

## See also

- [Concepts: Boundary conditions](../concepts/boundary_conditions.md)
- [Reference: Boundary](../reference/boundary.md)
