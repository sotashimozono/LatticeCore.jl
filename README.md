# LatticeCore.jl

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono.com/LatticeCore.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://codes.sota-shimozono.com/LatticeCore.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-v1.10+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

[![codecov](https://codecov.io/gh/sotashimozono/LatticeCore.jl/graph/badge.svg?token=oSYKPUteiH)](https://codecov.io/gh/sotashimozono/LatticeCore.jl)
[![Build Status](https://github.com/sotashimozono/LatticeCore.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sotashimozono/LatticeCore.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**LatticeCore** is the abstract interface layer for a family of
Julia packages that simulate physical systems on geometric lattices.
It defines the types and trait vocabulary that every lattice —
periodic or aperiodic, finite or conceptually infinite — must
implement, then ships a pair of reference implementations
(`LineLattice`, `SimpleSquareLattice`) and a momentum-space layer
so the interface is verifiable end-to-end without any other
dependency.

## Features

- `AbstractLattice{D, T}` with trait-based extension (`TopologyTrait`,
  `Periodic` / `Aperiodic`, `is_bipartite`, `reciprocal_support`,
  `size_trait`).
- **Per-axis boundary conditions** (`PeriodicAxis`, `OpenAxis`,
  `TwistedAxis`) composed into a `LatticeBoundary`, with a
  modifier slot for non-topological reweightings like SSD. Mixed
  BCs (cylinders) are first-class.
- **Coordinate system and indexing split**: `RealSpace`,
  `LatticeCoord`, `HigherDimCoord`, plus `RowMajor` / `ColMajor` /
  `Snake` indexing strategies decoupled from the coordinates they
  linearise.
- **Site types** (`IsingSite`, `PottsSite{Q}`, `XYSite`,
  `HeisenbergSite`, `EmptySite`) stored through three layouts
  (`UniformLayout`, `SublatticeLayout`, `ExplicitLayout`) so
  mixed-spin and disordered models compose cleanly. Geometric
  sublattice and physical site type live on separate axes on
  purpose.
- **Element centering** trait (`VertexCenter` by default,
  `BondCenter` / `PlaquetteCenter` / `CellCenter` as extension
  points) so dimer, gauge-link, and flux variables can be modelled
  without a breaking rewrite.
- **Momentum-space layer**: `AbstractMomentumLattice`,
  `PeriodicMomentumLattice`, `monkhorst_pack` / `gamma_centered`
  mesh constructors, and a naive `structure_factor`. Type
  skeletons (`HyperReciprocalLattice`, `BraggPeakSet`,
  `AcceptanceWindow`) for downstream quasicrystal Fourier
  analysis.
- **Lazy / infinite lattice hooks**: `materialize(abstract; cutoff)`
  for the infinite-abstract → finite-materialisation pattern, and
  `require_finite(lat)` for Monte Carlo entry-point guards.
- **Reference lattices** (`LineLattice`, `SimpleSquareLattice`)
  with full boundary-condition, site-type, and reciprocal-lattice
  support, so the interface can be driven end-to-end inside
  LatticeCore itself.

## Installation

LatticeCore targets Julia 1.10 or later.

```julia
using Pkg
Pkg.add("LatticeCore")
```

Until the General registry entry lands, install directly from
GitHub:

```julia
Pkg.add(url = "https://github.com/sotashimozono/LatticeCore.jl.git")
```

## A minimal example

```julia
using LatticeCore
using StaticArrays

# 4x4 periodic square lattice with default Ising sites.
lat = SimpleSquareLattice(4, 4, PeriodicAxis())

num_sites(lat)           # 16
position(lat, 5)         # SVector(1.0, 2.0)
neighbors(lat, 5)        # [6, 4, 9, 1]
site_type(lat, 5)        # IsingSite{Int8}()

# Cylinder: periodic in x, open in y.
cyl = SimpleSquareLattice(4, 4,
    LatticeBoundary((PeriodicAxis(), OpenAxis())))
periodicity(cyl)         # Aperiodic()

# Reciprocal lattice and a naive structure factor.
ml = reciprocal_lattice(lat)
state = ones(Int8, num_sites(lat))                       # ferromagnet
structure_factor(lat, state, SVector(0.0, 0.0))          # ≈ 16.0

# Guard Monte Carlo entry points against non-finite lattices.
require_finite(lat)  # no-op for finite lattices
```

See the [full documentation](https://codes.sota-shimozono.com/LatticeCore.jl/stable/)
for the concept columns (Bravais lattices, reciprocal space,
boundary conditions, quasiperiodic order, site types), the API
guide, and the auto-generated reference.

## Where LatticeCore fits

```text
LatticeCore.jl            (this package, abstract interface)
      │
      ├── Lattice2D.jl        (periodic 2D catalogue — downstream)
      ├── QuasiCrystal.jl     (cut-and-project quasicrystals — downstream)
      │
      └── Lattice2DMonteCarlo.jl  (classical MC runtime — downstream)
```

LatticeCore itself is deliberately small. It is the contract every
downstream package agrees on; it does not ship a production lattice
catalogue, a quasicrystal generator, or a Monte Carlo runtime.

## Quality assurance

The test suite runs [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl)
on every CI build (ambiguities, unbound type parameters, undefined
exports, stale deps, compat bounds, piracy). All checks pass on
`main`.

## License

MIT. See [`LICENSE`](LICENSE).
