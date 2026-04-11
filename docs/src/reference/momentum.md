# Momentum space API

See [Momentum space guide](../guide/momentum.md) for the narrative
version and
[Reciprocal lattice and Brillouin zone](../concepts/reciprocal_and_brillouin.md)
for the mathematics.

## Abstract types

```@docs
AbstractMomentumLattice
```

## Concrete momentum lattices

```@docs
PeriodicMomentumLattice
```

## Interface

```@docs
num_k_points
k_point
reciprocal_basis
```

## Mesh constructors

```@docs
monkhorst_pack
gamma_centered
```

## Entry points

```@docs
reciprocal_lattice
fourier_module
momentum_lattice
```

## Structure factor

```@docs
structure_factor
```

## Quasicrystal Fourier-module skeletons

```@docs
AcceptanceWindow
HyperReciprocalLattice
BraggPeakSet
```
