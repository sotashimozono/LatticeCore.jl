# Boundary API

See [Boundary conditions guide](../guide/boundary.md) and the
concept column [Boundary conditions](../concepts/boundary_conditions.md)
for the narrative background.

## Abstract types

```@docs
AbstractBoundaryCondition
AbstractAxisBC
AbstractBoundaryModifier
```

## Axis-level boundaries

```@docs
PeriodicAxis
OpenAxis
TwistedAxis
```

## Modifiers

```@docs
NoModifier
SSD
```

## Composite boundary

```@docs
LatticeBoundary
```

## Hook functions

```@docs
apply_axis_bc
axis_phase
bond_weight
```
