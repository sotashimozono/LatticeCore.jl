# Traits

LatticeCore uses small value types (Holy traits) to express optional
capabilities and invariants without inflating the concrete lattice's
type parameters. See [Lattice interface](../guide/lattice.md) for
where each trait is consumed.

## Topology

```@docs
TopologyTrait
topology
```

## Periodicity

```@docs
Periodic
Aperiodic
periodicity
is_bipartite
```

## Reciprocal-space support

```@docs
AbstractReciprocalSupport
HasReciprocal
HasFourierModule
NoReciprocal
reciprocal_support
```

## Size

```@docs
AbstractSizeTrait
FiniteSize
InfiniteSize
QuasiInfiniteSize
```
