# Lattices and unit cells

A **lattice** is the skeleton physicists use to describe a crystal:
an infinite, regular pattern of points in space obtained by
repeatedly translating a finite motif. Everything LatticeCore does —
storing positions, walking neighbours, computing structure factors —
is organised around this picture, so it is worth spending a few
minutes making the vocabulary precise.

## Bravais lattices

A **Bravais lattice** in ``d`` dimensions is the set of points

```math
\mathbf{R} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + \dots + n_d \mathbf{a}_d,
\qquad n_i \in \mathbb{Z},
```

where the **primitive vectors** ``\mathbf{a}_1, \dots, \mathbf{a}_d``
are linearly independent and fixed. The matrix whose columns are the
``\mathbf{a}_i`` is the **real-space basis** ``A``. In LatticeCore a
concrete lattice exposes it via `basis_vectors(lat)`.

Two Bravais lattices are the same if and only if they share a basis
up to an integer (unimodular) linear transformation. That means the
basis is never unique — the honest square lattice

```math
\mathbf{a}_1 = (1, 0), \quad \mathbf{a}_2 = (0, 1)
```

is the same Bravais lattice as the skewed choice

```math
\mathbf{a}_1' = (1, 0), \quad \mathbf{a}_2' = (1, 1),
```

and the code treats both as equivalent inputs. What *is* canonical
is the set of lattice points ``\{\mathbf{R}\}`` and everything that
can be computed from it: neighbour distances, reciprocal basis
(see [Reciprocal lattice and Brillouin zone](reciprocal_and_brillouin.md)),
density.

Two-dimensional Bravais lattices come in exactly five types
(oblique, rectangular, centred rectangular, hexagonal, square),
three-dimensional ones in fourteen. LatticeCore's reference
implementations ([`SimpleSquareLattice`](@ref),
[`LineLattice`](@ref)) are the square and 1D cases with unit
spacing, which is enough to exercise every piece of the interface.

## Unit cells and sublattices

A crystal is usually not *just* a Bravais lattice. Each lattice
point carries a **basis** — a finite cluster of atoms, or (for us) a
finite cluster of sites. Copies of that cluster sit at every lattice
point, and together they make the full crystal:

```math
\mathrm{Crystal} = \mathrm{Bravais\ lattice} + \mathrm{motif}.
```

The cluster is called the **unit cell**, and each member of the
cluster is a **sublattice**. A honeycomb lattice is a triangular
Bravais lattice with a two-site motif — one A site, one B site — and
that is exactly why the A and B sublattices never coincide.

LatticeCore records sublattice membership inside the
[`LatticeCoord`](@ref) type:

```julia
struct LatticeCoord{D}
    cell::NTuple{D, Int}   # which Bravais cell
    sublattice::Int        # which site inside the cell
end
```

A coordinate like `LatticeCoord((3, 2), 1)` means "the A site of the
(3, 2) cell". This lets the same cell index talk to multiple
sublattices without any tuple juggling.

### Sublattice vs site type

There are *two* useful ways to label a site, and they are
genuinely different:

- **Geometric sublattice** — *where* in the unit cell the site sits
  (A, B, C, ...). This is pure geometry and it is what
  `LatticeCoord.sublattice` stores.
- **Site type** — *what* physical degree of freedom lives on the
  site (Ising spin, XY rotor, Heisenberg vector, vacancy, gauge
  link, ...). This is captured by [`AbstractSiteType`](@ref) and
  [`AbstractSiteLayout`](@ref).

Conflating the two is the single most common modelling bug in a
mixed-spin simulation. LatticeCore keeps them on separate axes
deliberately: a honeycomb lattice can have an Ising site on the A
sublattice and an XY site on the B sublattice, and the code never
has to decide which field to ask. See
[Site types and spin models](site_types_and_spins.md) for the site
type side.

## Periodicity, bipartiteness, and coordination number

Three numbers describe a Bravais-like lattice at a glance:

- The **coordination number** ``z`` — how many nearest neighbours a
  site has. It is a topology invariant (square lattice: 4,
  triangular: 6, honeycomb: 3, kagome: 4), visible in LatticeCore
  through `length(neighbors(lat, i))` for an interior site.
- The **bipartiteness** of the neighbour graph. A lattice is
  bipartite if the sites can be two-coloured so that every bond
  connects different colours. Square and honeycomb are bipartite;
  triangular is not. Bipartiteness controls whether a Neel
  antiferromagnetic ground state even *exists* and is exposed as
  [`is_bipartite`](@ref).
- The **periodicity** of the underlying crystal. A Bravais crystal
  is periodic by construction; a Penrose tiling is *aperiodic* even
  though it is highly ordered. LatticeCore captures this with the
  [`Periodic`](@ref) / [`Aperiodic`](@ref) trait attached to every
  lattice.

## Why this matters for the API

Once these concepts are named, the LatticeCore interface almost
writes itself:

- `position(lat, i)` returns a real-space point — a sum of
  primitive vectors plus a sublattice offset. It is the one function
  that carries the full geometric information of the lattice.
- `neighbors(lat, i)` encodes the topology — which *other* sites
  are connected by a bond, independent of positions.
- `boundary(lat)` decides what happens at the edges of a **finite
  piece** of the otherwise infinite Bravais lattice. That is the
  subject of [Boundary conditions](boundary_conditions.md).

The rest of the package (site types, reciprocal lattice, structure
factor, lazy materialisation) is built out of these three calls.

## Further reading

- N. W. Ashcroft and N. D. Mermin, *Solid State Physics*, Ch. 4–5
  (Bravais lattices and reciprocal lattices).
- C. Kittel, *Introduction to Solid State Physics*, Ch. 1–2.
- S. Sachdev, *Quantum Phase Transitions*, Ch. 1 (for the role of
  lattice geometry in ordering transitions).
