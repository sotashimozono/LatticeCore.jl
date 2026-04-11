# Coordinate systems and indexing

LatticeCore distinguishes two related-but-independent concerns:

1. **Coordinate systems** — how a site's *position* is described.
   A point has three natural descriptions: Cartesian real-space,
   lattice-cell + sublattice, and higher-dimensional hyper
   coordinates (for cut-and-project quasicrystals).
2. **Indexing** — how a coordinate is linearised into a 1-based
   integer site index. This is the rule that says "the A site of
   cell (2, 3) is site number 7".

Keeping them separate means a lattice can switch its indexing
without touching its coordinate representation, and vice versa.

## The three coordinate kinds

```julia
abstract type AbstractCoordinate{D} end
```

Every concrete subtype describes a point in the same ``D``-
dimensional physical space — just in a different chart.

### RealSpace

Cartesian coordinates, backed by `StaticArrays.SVector`:

```julia
rs = RealSpace((1.5, 2.0))          # from a Tuple
rs = RealSpace(SVector(1.5, 2.0))   # from an SVector
rs.x                                 # SVector{2, Float64}([1.5, 2.0])
```

`RealSpace{D, T}` is parameterised on the storage type so you can
choose `Float32` for memory pressure or `BigFloat` for precision.

### LatticeCoord

Integer cell index plus a 1-based **geometric** sublattice id:

```julia
lc = LatticeCoord((3, 4))           # sublattice defaults to 1
lc = LatticeCoord((3, 4), 2)        # explicit sublattice (B site)
lc.cell         # (3, 4)
lc.sublattice   # 2
```

The sublattice field is **geometric** — where in the unit cell the
site sits — not physical. A honeycomb A/B distinction uses this
field; a mixed-spin Ising-on-A / XY-on-B model separately uses
[`SublatticeLayout`](@ref) to attach site types to these sublattices.
See [Site types and layouts](site_type.md).

### HigherDimCoord

For cut-and-project quasicrystals, the natural description of a
point is as a lattice point in a *higher*-dimensional host
lattice:

```julia
h = HigherDimCoord{2}(SVector(1.0, 2.0, 3.0, 4.0, 5.0))   # Penrose: 5D → 2D
h isa HigherDimCoord{2, 5, Float64}   # true
```

`DPhys` (the first type parameter) is the physical dimension;
`DHyper` is the host-lattice dimension. The projection from
`DHyper` to `DPhys` is the lattice's responsibility — it lives in
the quasicrystal implementation, not in the coordinate type.

## Conversions

All three kinds talk to each other through three generic
functions, each dispatched on the concrete lattice type:

```julia
to_real(lat, coord)     # → RealSpace
to_lattice(lat, coord)  # → LatticeCoord
to_hyper(lat, coord)    # → HigherDimCoord
```

Example on the reference square lattice:

```julia
lat = SimpleSquareLattice(3, 3)

rs = to_real(lat, LatticeCoord((2, 3)))
rs.x                                    # SVector(2.0, 3.0)

lc = to_lattice(lat, RealSpace((2.0, 3.0)))
lc.cell                                 # (2, 3)
lc.sublattice                           # 1

# Round-trip
to_lattice(lat, to_real(lat, LatticeCoord((2, 3))))
# → LatticeCoord((2, 3), 1)
```

Identity conversions are defined for every coordinate kind, so
you can call `to_real(lat, coord)` blindly even if `coord` is
already a `RealSpace` — it will be returned unchanged.

For lattices that do not live in a higher-dimensional host (every
periodic lattice), `to_hyper` simply is not defined and calling it
will raise `MethodError`. The reciprocal is also true: a
quasicrystal may not implement `to_lattice(lat, ::RealSpace)` if
the map is not unique.

## Indexing strategies

LatticeCore ships three indexings out of the box, all under
[`AbstractIndexing`](@ref):

- [`RowMajor`](@ref) — C-style. `x` is the fast axis within a row.
- [`ColMajor`](@ref) — Fortran-style. `y` is the fast axis.
- [`Snake`](@ref) — row-major but with alternating row direction.

Each ships with `site_index` and `lattice_coord` methods for 1D
and 2D. The sublattice is the **innermost** index, so `nsub = 1`
reduces to the usual single-sublattice formulas and can be read
at a glance.

### RowMajor on a 3 × 2 lattice

```julia
# Layout:
#   row y=1:  sites 1 2 3
#   row y=2:  sites 4 5 6
site_index(RowMajor(), (3, 2), 1, LatticeCoord((1, 1)))   # 1
site_index(RowMajor(), (3, 2), 1, LatticeCoord((3, 1)))   # 3
site_index(RowMajor(), (3, 2), 1, LatticeCoord((1, 2)))   # 4
site_index(RowMajor(), (3, 2), 1, LatticeCoord((3, 2)))   # 6

lattice_coord(RowMajor(), (3, 2), 1, 4).cell   # (1, 2)
```

### ColMajor on the same dims

```julia
# Layout:
#   col x=1:  sites 1 2
#   col x=2:  sites 3 4
#   col x=3:  sites 5 6
site_index(ColMajor(), (3, 2), 1, LatticeCoord((1, 1)))   # 1
site_index(ColMajor(), (3, 2), 1, LatticeCoord((1, 2)))   # 2
site_index(ColMajor(), (3, 2), 1, LatticeCoord((2, 1)))   # 3
```

### Snake

Row-major with every second row walked backwards. For a 3 × 3
snake, the cell index inside each row is flipped when `y` is
even:

```text
y = 1 (forward):   1 2 3
y = 2 (reverse):   6 5 4
y = 3 (forward):   7 8 9
```

```julia
site_index(Snake(), (3, 3), 1, LatticeCoord((1, 1)))   # 1
site_index(Snake(), (3, 3), 1, LatticeCoord((3, 2)))   # 4
site_index(Snake(), (3, 3), 1, LatticeCoord((1, 2)))   # 6
```

Snake is useful when physical adjacency between neighbouring
indices matters — serialisation layouts, cache-friendly sweeps,
some SIMD-friendly scheduling tricks.

## Round-trip correctness

`site_index` and `lattice_coord` are **exact inverses** for
every indexing: the LatticeCore test suite runs an exhaustive
check for every `(indexing, dims, nsub)` combination, hitting
every site index exactly once. That is the contract for every new
indexing you might write.

## Adding a custom indexing

Adding a new indexing strategy is two method definitions:

```julia
struct HilbertOrder <: LatticeCore.AbstractIndexing end

function LatticeCore.site_index(
    ::HilbertOrder, dims::NTuple{2, Int}, nsub::Int, coord::LatticeCoord{2}
)
    # ... map (cx, cy) into Hilbert curve coordinate, then a 1-based index ...
end

function LatticeCore.lattice_coord(
    ::HilbertOrder, dims::NTuple{2, Int}, nsub::Int, i::Int
)
    # ... inverse Hilbert curve map ...
end
```

Once those are in place, any lattice that routes its cell
bookkeeping through `site_index` / `lattice_coord` can switch
indexing strategies by instance.

## See also

- [Concepts: Lattices and unit cells](../concepts/lattice_and_unit_cells.md)
- [Reference: Coordinates](../reference/coordinates.md)
- [Reference: Indexing](../reference/indexing.md)
