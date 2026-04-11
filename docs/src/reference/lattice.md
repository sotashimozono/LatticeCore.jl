# AbstractLattice API

Auto-generated reference for the `AbstractLattice` interface. See
[Lattice interface](../guide/lattice.md) for the narrative version
and [Lattices and unit cells](../concepts/lattice_and_unit_cells.md)
for the physical background.

## Abstract type

```@docs
AbstractLattice
```

## Required interface

```@docs
num_sites
position
neighbors
boundary
size_trait
site_layout
```

## Derived / default helpers

```@docs
dimension
positions
is_finite
num_sublattices
sublattice
```
