# Boundary conditions

LatticeCore describes a boundary condition as a composite of two
things: a tuple of **axis boundary conditions** (one per spatial
axis) and a single non-topological **modifier**. Axis BCs decide
which candidate bonds exist; modifiers decide how existing bonds
are weighted.

For the physics behind the taxonomy, see the concept column
[Boundary conditions](../concepts/boundary_conditions.md). This
guide is the "how to use it" side.

## The composite

```julia
LatticeBoundary(
    axes     = (PeriodicAxis(), OpenAxis()),
    modifier = NoModifier(),
)
```

- `axes::NTuple{D, <:AbstractAxisBC}` — exactly one per physical
  axis. Legal values: [`PeriodicAxis`](@ref), [`OpenAxis`](@ref),
  [`TwistedAxis(θ)`](@ref TwistedAxis).
- `modifier::AbstractBoundaryModifier` — defaults to
  [`NoModifier`](@ref); [`SSD`](@ref) is the other type LatticeCore
  currently ships.

A uniform choice is usually easier to read:

```julia
LatticeBoundary((PeriodicAxis(), PeriodicAxis()))    # 2D torus
LatticeBoundary((OpenAxis(), OpenAxis()))            # 2D open sample
```

and the reference lattice constructors even accept a single
`AbstractAxisBC` as shorthand, which they broadcast to every axis:

```julia
SimpleSquareLattice(3, 3, PeriodicAxis())   # 3×3 torus
SimpleSquareLattice(3, 3, OpenAxis())       # 3×3 open sample
```

Mixed BCs have to go through `LatticeBoundary` explicitly.

## Worked example: a 3×4 cylinder

A cylinder is periodic along one axis and open along the other.
Here we make a 3×4 square cylinder with PBC along `x` and OBC
along `y`, and we read the neighbours of every site:

```julia
using LatticeCore

cyl = SimpleSquareLattice(3, 4,
    LatticeBoundary((PeriodicAxis(), OpenAxis())))

for i in 1:num_sites(cyl)
    println(i, "  →  ", neighbors(cyl, i))
end
```

What comes out:

```text
1  →  [2, 3, 4]          # corner: +x, -x wraps to site 3, +y
2  →  [3, 1, 5]
3  →  [1, 2, 6]          # right edge: +x wraps to site 1
4  →  [5, 6, 7, 1]       # interior x-edge, still boundary in y
5  →  [6, 4, 8, 2]
6  →  [4, 5, 9, 3]
7  →  [8, 9, 10, 4]
8  →  [9, 7, 11, 5]
9  →  [7, 8, 12, 6]
10 →  [11, 12, 7]        # top-row corner: no +y neighbour
11 →  [12, 10, 8]
12 →  [10, 11, 9]
```

Along ``x`` every site has exactly two neighbours (the wrap
closes it into a 3-cycle). Along ``y`` only the interior rows
have both neighbours; the top and bottom rows drop the
out-of-range neighbour. This is what a cylinder looks like from
the neighbour list.

[`periodicity`](@ref) reports `Aperiodic()` for the same lattice
because the overall structure is not fully translation-invariant,
even though one axis wraps. If you want a predicate on a single
axis, look directly at `cyl.boundary.axes[i]`.

## Checking wrapped bond vectors

PBC is subtle when you care about the displacement vector of a
boundary-crossing bond. The reference lattices override
[`neighbor_bonds`](@ref) so that a wrapped bond carries the
**shortest** displacement (``\pm 1`` in each axis), not the raw
``\mathrm{position}(j) - \mathrm{position}(i)``:

```julia
lat = SimpleSquareLattice(3, 3, PeriodicAxis())
nb1 = collect(neighbor_bonds(lat, 1))

wrap_x = first(b for b in nb1 if b.j == 3)      # wraps in x
wrap_x.vector                                    # SVector(-1.0, 0.0)

wrap_y = first(b for b in nb1 if b.j == 7)      # wraps in y
wrap_y.vector                                    # SVector(0.0, -1.0)
```

Without the override the displacements would be `(2, 0)` and
`(0, 2)`, which breaks bond-centred observables and SSD weights.
If you write a custom lattice, overriding `neighbor_bonds` to emit
the wrapped vector is the canonical thing to do.

## The hook layer

Every axis BC is a dispatch point for two generic functions:

- [`apply_axis_bc(axis_bc, idx, L)`](@ref apply_axis_bc) —
  returns `(wrapped, is_valid)`. The connectivity hook.
- [`axis_phase(axis_bc, idx, L)`](@ref axis_phase) — returns a
  `ComplexF64`. Only [`TwistedAxis`](@ref) returns a non-unit
  phase; the classical MC code paths are free to ignore it.

Modifiers dispatch through [`bond_weight(modifier, lat, i, j)`](@ref
bond_weight). [`NoModifier`](@ref) returns `1.0` unconditionally;
[`SSD`](@ref)'s weighting function is a stub until the MC layer
lands a `center(lat)` accessor.

## Adding a new axis BC

Adding a BC kind is three short methods — you do not edit any
existing LatticeCore file. For instance, an **anti-periodic**
boundary (``\theta = \pi``) is just

```julia
struct AntiPeriodicAxis <: LatticeCore.AbstractAxisBC end

LatticeCore.apply_axis_bc(::AntiPeriodicAxis, idx::Int, L::Int) =
    (mod1(idx, L), true)

LatticeCore.axis_phase(::AntiPeriodicAxis, idx::Int, L::Int) =
    (idx < 1 || idx > L) ? complex(-1.0, 0.0) : complex(1.0, 0.0)
```

You can then feed it straight into `LatticeBoundary`:

```julia
LatticeBoundary((AntiPeriodicAxis(), PeriodicAxis()))
```

and every concrete lattice that walks its bonds through
`apply_axis_bc` — including the reference square lattice — picks
it up without any further code change. This is the extension
mechanism the architecture assumes: one new type, one or two
method definitions, no core edits.

## See also

- [Concepts: Boundary conditions](../concepts/boundary_conditions.md)
- [Reference: Boundary](../reference/boundary.md)
