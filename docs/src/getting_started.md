# Getting Started

## Installation

LatticeCore is not yet registered in the Julia General registry.
Install it directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/sotashimozono/LatticeCore.jl.git")
```

## A minimal example

Every downstream package — periodic lattices, quasicrystals, the
Monte Carlo runtime — builds on top of `AbstractLattice`. LatticeCore
itself ships two concrete reference implementations so you can
exercise the full stack without any other dependency.

```julia
using LatticeCore
using StaticArrays

# 1D periodic chain of 8 sites with default Ising site type
chain = LineLattice(8, PeriodicAxis())

num_sites(chain)            # 8
position(chain, 1)          # SVector(1.0)
neighbors(chain, 1)         # [8, 2]       (wraps to the far end)
site_type(chain, 1)         # IsingSite{Int8}()

# 4x4 periodic square lattice
sq = SimpleSquareLattice(4, 4, PeriodicAxis())

num_sites(sq)               # 16
position(sq, 5)             # SVector(1.0, 2.0)
neighbors(sq, 5)             # [6, 4, 9, 1]
is_bipartite(sq)            # true (both axis lengths are even)

# Cylinder: periodic in x, open in y
cyl = SimpleSquareLattice(4, 4,
    LatticeBoundary((PeriodicAxis(), OpenAxis())))
periodicity(cyl)            # Aperiodic()  (one axis is open)

# Reciprocal lattice and a naive structure factor
ml = reciprocal_lattice(sq)
state = Int8[1 for _ in 1:num_sites(sq)]           # ferromagnet
structure_factor(sq, state, SVector(0.0, 0.0))     # ≈ 16.0
```

## The concepts

A lattice in LatticeCore is described along five independent axes:

1. **Topology** — which candidate bonds exist between sites.
   Captured by the topology trait and the lattice's native
   neighbour function.
2. **Boundary** — how those candidate bonds wrap (or don't) at the
   edges of the lattice. Expressed as a
   [`LatticeBoundary`](@ref) composed of per-axis
   [`AbstractAxisBC`](@ref) and an optional
   [`AbstractBoundaryModifier`](@ref).
3. **Coordinate system** — how a physical position is identified:
   as a [`RealSpace`](@ref) point, a [`LatticeCoord`](@ref) (cell
   index + geometric sublattice), or a [`HigherDimCoord`](@ref) for
   cut-and-project quasicrystals.
4. **Geometric sublattice** — the A/B/C role of a site within its
   unit cell. Lives in `LatticeCoord.sublattice`.
5. **Site type** — the physical degree of freedom on the site
   ([`IsingSite`](@ref), [`XYSite`](@ref), [`HeisenbergSite`](@ref),
   ...), stored via an [`AbstractSiteLayout`](@ref)
   ([`UniformLayout`](@ref), [`SublatticeLayout`](@ref),
   [`ExplicitLayout`](@ref)).

These are deliberately **orthogonal**: you can put any site type on
any layout on any boundary on any topology, and they compose via
multiple dispatch without changing core types.

## Guarding Monte Carlo entry points

Monte Carlo algorithms can only run on finite lattices. Guard the
entry point of your `run!` function with [`require_finite`](@ref):

```julia
function run!(rng, state, lat::AbstractLattice, model, alg; kwargs...)
    require_finite(lat)
    # ... safe to iterate 1:num_sites(lat) from here ...
end
```

For a conceptually-infinite structure, lower it to a finite lattice
first with [`materialize`](@ref):

```julia
abstract_fibonacci = InfiniteFibonacci(...)          # downstream package
lat = materialize(abstract_fibonacci; depth = 12)    # FiniteSize
run!(rng, state, lat, model, alg)
```

## Where to go next

- [AbstractLattice interface](guide/lattice.md)
- [Boundary conditions](guide/boundary.md)
- [Coordinate systems](guide/coordinates.md)
- [Site types and layouts](guide/site_type.md)
- [Momentum space](guide/momentum.md)
- [Lazy / infinite lattices](guide/lazy_infinite.md)
