# Translation-cell / orbit API

The translation-cell layer describes a lattice by its unit-cell motif —
a finite site basis and a finite bond motif — plus a translation action,
and exposes the fundamental domain through *orbit* accessors that work
for finite and infinite lattices alike. See the
[Lazy / infinite lattices guide](../guide/lazy_infinite.md) for context.

## Motif identities

```@docs
CellSite
CellBond
```

## Motif interface

```@docs
translation_vectors
num_basis_sites
basis_position
cell_bonds
```

## Orbit accessors

```@docs
site_orbits
bond_orbits
plaquette_orbits
element_orbits
```

## Lazy access and orbit positions

```@docs
cell_position
incident_cell_bonds
neighbors_at
element_orbit_position
```
