# Site types and layouts

A **site type** is a value-level description of the physical
degree of freedom living on a lattice site. A **site layout** is
how those site types are stored per-lattice. LatticeCore keeps the
two concerns distinct from the *geometric* sublattice
(`LatticeCoord.sublattice`) so mixed-spin and diluted models
compose naturally.

For the physics behind the taxonomy, see the concept column
[Site types and spin models](../concepts/site_types_and_spins.md).

## The built-in site types

All of the classical spin flavours are available out of the box:

| Site type                | State                    | Default storage     |
| ------------------------ | ------------------------ | ------------------- |
| [`IsingSite`](@ref)      | ``\pm 1``                | `Int8`              |
| [`PottsSite{Q}`](@ref)   | `1..Q`                   | `Int8`              |
| [`XYSite`](@ref)         | angle in `[0, 2π)`       | `Float64`           |
| [`HeisenbergSite`](@ref) | unit vector on ``S^2``   | `SVector{3, Float64}` |
| [`EmptySite`](@ref)      | `Nothing`                | `Nothing`           |

Each one is constructed by calling the struct with no arguments
for the default storage type:

```julia
IsingSite()           # IsingSite{Int8}()
IsingSite{Int}()      # Int storage (more memory, wider range)

PottsSite(3)          # 3-state Potts, Int8 storage
XYSite()              # Float64 angle
HeisenbergSite()      # SVector{3, Float64} unit vector
EmptySite()           # vacancy
```

The required interface on every site type is

```julia
state_type(st)                  # the Julia type storing one state
random_state(rng, st)           # uniform sampler for initialisation
```

plus the optional

```julia
zero_state(st)      # canonical zero state, if meaningful
domain(st)          # finite iterator over the state space, for discrete types
```

and the element-centering trait

```julia
element_type(st)    # default VertexCenter()
```

which is covered in more detail below.

## Picking a layout

Sites with site types are stored in the lattice's
[`AbstractSiteLayout`](@ref). Three concrete layouts ship:

### UniformLayout — every site the same

```julia
sq = SimpleSquareLattice(4, 4)                 # defaults to IsingSite
site_layout(sq)                                # UniformLayout{IsingSite{Int8}}
site_type(sq, 1) === site_type(sq, 16)         # true

# Swap in a different site type via the keyword argument
heis = SimpleSquareLattice(4, 4; layout = UniformLayout(HeisenbergSite()))
site_type(heis, 1)                              # HeisenbergSite{Float64}()
```

`UniformLayout` has **zero per-site memory overhead** — the single
site type is stored once on the lattice, so the MC hot loop can
treat `site_type(lat, i)` as a loop-invariant and constant-fold it.
This is the layout you want for homogeneous Ising / XY / Heisenberg
models, i.e. the overwhelming majority.

### SublatticeLayout — one site type per geometric sublattice

```julia
# 6 sites with the A/B/A/B/A/B pattern of a chain
mixed = SublatticeLayout(
    (IsingSite(), XYSite()),      # per-sublattice site types
    [1, 2, 1, 2, 1, 2],           # sublattice id of each site
)

site_type(mixed, 1) isa IsingSite   # true
site_type(mixed, 2) isa XYSite      # true
site_type(mixed, 5) isa IsingSite   # true
```

The first field is an `NTuple` of site type instances, one per
geometric sublattice. The second is a `Vector{Int}` mapping every
site index to its sublattice id. Downstream code can then write
polymorphic interactions:

```julia
interaction(m::MixedSpinModel, ::IsingSite, ::XYSite, bond, s, θ) =
    -m.J_AB * s * cos(θ)
```

and the multiple-dispatch mechanism picks out every A–B bond
automatically.

### ExplicitLayout — quenched per-site heterogeneity

```julia
types = AbstractSiteType[IsingSite(), IsingSite(), XYSite(), HeisenbergSite()]
disorder = ExplicitLayout(types)

site_type(disorder, 1) isa IsingSite       # true
site_type(disorder, 3) isa XYSite          # true
site_type(disorder, 4) isa HeisenbergSite  # true
```

`ExplicitLayout` stores one site type per site, which is the right
representation for quenched disorder, defect models, or any
configuration where the pattern cannot factor through a small
sublattice tuple. The per-site storage cost is the trade-off.

## A custom site type from scratch

Custom site types are what you add when the existing flavours do
not fit — for instance, a Blume–Capel spin (``S = 1``), a
``U(1)`` clock rotor, a Kitaev bond variable, or a dimer variable.
The required steps are

1. Define a subtype of [`AbstractSiteType`](@ref).
2. Implement `state_type` and `random_state`.
3. (Optional) implement `zero_state`, `domain`, `element_type`.

A three-state Blume–Capel example:

```julia
using Random

struct BlumeCapelSite <: LatticeCore.AbstractSiteType end

LatticeCore.state_type(::BlumeCapelSite) = Int8
LatticeCore.random_state(rng, ::BlumeCapelSite) =
    Int8(rand(rng, (-1, 0, 1)))
LatticeCore.zero_state(::BlumeCapelSite) = Int8(0)
LatticeCore.domain(::BlumeCapelSite) = (Int8(-1), Int8(0), Int8(1))
```

It plugs straight into the existing layouts:

```julia
lat = SimpleSquareLattice(4, 4; layout = UniformLayout(BlumeCapelSite()))
site_type(lat, 1) isa BlumeCapelSite     # true
```

## Element centering

Every `AbstractSiteType` carries an `element_type` trait declaring
where its DOF lives. The default is [`VertexCenter`](@ref), but
bond / plaquette / cell-centered variables are first-class through
this trait:

```julia
struct DimerVariable <: LatticeCore.AbstractSiteType end
LatticeCore.element_type(::DimerVariable) = BondCenter()
LatticeCore.state_type(::DimerVariable) = Bool
LatticeCore.random_state(rng, ::DimerVariable) = rand(rng, Bool)
```

A downstream Monte Carlo layer that wants to know whether a site
type lives on a vertex or a bond asks `element_type(st)`. For the
homogeneous classical spin models that currently dominate the
LatticeCore test suite, this trait returns `VertexCenter()` and
every code path ignores it, so adding a bond-centered type does
not bloat vertex-only code.

See [Concepts: Quasiperiodic order](../concepts/quasiperiodic_order.md)
for why this matters when you start thinking about dimer models,
gauge theories, and flux variables.

## Querying a lattice

Three accessors on `AbstractLattice` let downstream code read the
site structure of a lattice:

```julia
site_layout(lat)              # AbstractSiteLayout
site_type(lat, i)             # AbstractSiteType at site i
num_sublattices(lat)          # Int, defaults to 1
sublattice(lat, i)            # geometric sublattice id, defaults to 1
```

The default implementations of `site_type(lat, i)` delegate to
`site_layout(lat)`, so a concrete lattice only needs to implement
`site_layout` to get everything for free. `num_sublattices` and
`sublattice` default to `1` — lattices with non-trivial geometric
sublattices (honeycomb, kagome, Lieb) override them.

## See also

- [Concepts: Site types and spin models](../concepts/site_types_and_spins.md)
- [Reference: Site type](../reference/site_type.md)
- [Reference: Site layout](../reference/site_layout.md)
- [Reference: Lattice element](../reference/lattice_element.md)
