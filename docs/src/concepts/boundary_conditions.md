# Boundary conditions

Every lattice simulation you will ever run is finite. The real
crystal you are approximating, on the other hand, is almost always
huge compared to your sample — so the *boundary* of your sample is
where the approximation bleeds through. Picking a boundary condition
is picking the way that bleed is controlled, and it changes both the
physics you measure and the numerics you use to measure it.

LatticeCore's answer is to give you enough vocabulary to pick the
right one on purpose, and to compose it axis by axis.

## The four regimes

Four boundary conditions appear over and over again in lattice
simulations. They differ in how the lattice graph is closed off
(topology) and in whether the bonds that close it carry a phase or
a weight.

### Periodic (torus)

Identify each axis with its opposite edge. The sample becomes a
torus, the neighbour graph is translation-invariant, and every site
looks locally identical to a bulk site. This is the gentlest finite-
size approximation for ground-state physics:

- Translation symmetry is preserved, so crystal momentum is a good
  quantum number and the **discrete reciprocal lattice** is exact.
- Finite-size corrections decay like ``1/L^d`` or better for
  gapped systems — the Polyakov loops of the problem still wrap,
  but their effect is uniform across the sample.
- Monte Carlo algorithms integrate the ergodic sector directly, no
  edge carve-outs required.

In LatticeCore, `PeriodicAxis()` on every axis plus `NoModifier()`
is what `SimpleSquareLattice(Lx, Ly)` gives you out of the box.

### Open (real surface)

Don't close the axis at all. Edge and corner sites have fewer
neighbours than bulk sites, the translation symmetry is broken, and
the system *has* a genuine surface — with all the surface physics
that entails:

- Open sample calculations carry surface free energies and surface
  modes that are physically real and often interesting in their
  own right (surface plasmons, edge states, Kondo surface effects).
- Friedel oscillations near the boundary mean that "bulk" averages
  should be taken over an interior window, not the whole sample.
- DMRG / tensor-network algorithms converge much faster on OBC
  than on PBC in 1D because the ends of the chain have low
  entanglement.

`OpenAxis()` on every axis switches the reference lattices into
this regime.

### Cylinders and strips

Mix the two: periodic along one axis, open along another. On a
square you get a cylinder, on a cubic system you get a tube. This
is the workhorse geometry for tensor networks because it combines
the small edge effects of PBC along one direction with the clean
finite-width convergence of OBC along the orthogonal one. It is
also the natural choice for studying 1D-like quasi-long-range
order on top of a 2D lattice:

- Momentum along the periodic axis remains a good quantum number;
  momentum along the open axis does not.
- The circumference `Lx` sets the IR cutoff along the periodic axis
  (gaps close like `1/Lx` in gapless phases).
- On a honeycomb cylinder, the choice of circumferential axis
  (zigzag vs armchair) changes the edge physics qualitatively.

LatticeCore's `LatticeBoundary((PeriodicAxis(), OpenAxis()))` on a
`SimpleSquareLattice` is an explicit cylinder.

### Twisted periodic (flux and Bloch phases)

Close the sample like a torus, but attach a phase ``e^{i\theta}``
to any bond that crosses the boundary. Physically, this is **what a
lattice model sees of a magnetic flux**: the Peierls substitution
that takes a hopping integral ``t_{ij}`` to
``t_{ij} e^{i \oint \mathbf{A} \cdot d\mathbf{\ell}}`` localises the
flux on the boundary-crossing bonds of a sufficiently clever gauge
choice. The lattice topology is identical to PBC — every site has
the same number of neighbours — but the Hamiltonian is different.

Twisted BC are how lattice simulations talk about:

- Magnetic flux quantisation (Byers–Yang). Threading a flux
  ``\theta = \pi/2, \pi, 3\pi/2`` and tracking the ground-state
  energy picks out a superfluid response or a topological
  invariant.
- Finite-size extrapolation to zero-temperature conductance in
  tight-binding chains.
- Averaging over Bloch phases as a stand-in for momentum integration.

`TwistedAxis(θ)` stores the twist angle. The connectivity function
[`apply_axis_bc`](@ref) treats it identically to `PeriodicAxis`
(the neighbour graph is the same), while [`axis_phase`](@ref)
returns the phase factor that bond evaluation picks up.

## Why per-axis?

Many physics codes expose boundary conditions as a single global
flag: "PBC" or "OBC". That design breaks as soon as you want a
cylinder or a twisted flux — the axes stop being symmetric, and
you have to introduce per-axis escape hatches anyway.

LatticeCore follows the convention used by DMRG and tensor-network
libraries: store **one [`AbstractAxisBC`](@ref) per axis**, plus an
optional [`AbstractBoundaryModifier`](@ref) that reweights bonds
non-topologically, in a composite [`LatticeBoundary`](@ref). The
reference square lattice dispatches through each axis independently:

```julia
function neighbors(l::SimpleSquareLattice, i::Int)
    x, y = _site_to_xy(l, i)
    bx, by = l.boundary.axes

    ns = Int[]
    # +x direction
    nx, ok = apply_axis_bc(bx, x + 1, l.Lx)
    ok && push!(ns, _xy_to_site(l, nx, y))
    # ... and so on for −x, +y, −y
    return ns
end
```

Because `apply_axis_bc` dispatches on the axis BC's concrete type,
the compiler strips away the wrong-BC branches and the code is
specialised for the specific combination you actually use.

## Modifiers, and why SSD lives outside the axis tuple

Topological boundary conditions (the axis tuple) decide which
bonds *exist*. **Modifiers** are the other story: a non-topological
reweighting applied on top of an existing bond graph.

The archetypal modifier is **sine-square deformation** (SSD),
introduced in Gendiar, Krcmar & Nishino (2009). It multiplies
every bond's contribution to the Hamiltonian by the envelope

```math
f(\mathbf{r}) = \sin^2\!\left(\pi \, \frac{|\mathbf{r} - \mathbf{r}_{\text{centre}}|}{L}\right),
```

which is zero at the sample edges and maximal in the bulk. Open
boundary conditions dressed with an SSD envelope reproduce the
ground-state correlation functions of *periodic* systems almost
perfectly, because the smooth envelope suppresses surface
excitations before they can contaminate the bulk. This is
surprising and extremely useful; it is also obviously not a
wrapping rule.

LatticeCore's [`NoModifier`](@ref) and [`SSD`](@ref) live in a
separate slot of [`LatticeBoundary`](@ref). You can combine any
axis-level BC with any modifier without touching the core code:

```julia
LatticeBoundary((PeriodicAxis(), OpenAxis()), SSD(10.0))
```

[`bond_weight(::SSD, lat, i, j)`](@ref bond_weight) evaluates the
multi-axis envelope `∏_d sin²(π (c_d − 1/2) / L_d)` at each endpoint
and averages, reading the per-axis lengths `L_d` from the lattice's
[`size_trait`](@ref). The `L` field on `SSD(L)` is informational and
not used by the canonical evaluation; downstream packages may
override the method on a more specific lattice type to supply
alternative scales.

## Finite-size scaling and BC choice

A boundary condition is not just a technicality — it changes the
leading finite-size corrections to everything you measure. Three
rules of thumb:

1. **Gapped phases** have exponentially small PBC corrections and
   algebraic OBC corrections. Use PBC if you want to extract the
   bulk gap from a small sample.
2. **Critical points** have universal finite-size scaling. Usually
   you want PBC because it preserves translation symmetry, but
   cylinders are standard in 2D because they give you a natural
   bond dimension cutoff for DMRG.
3. **Surface physics** — edge modes, surface magnetism, surface
   plasmons — obviously needs OBC. Cylinders let you isolate one
   direction of surface physics at a time.

The [`periodicity`](@ref) trait on a LatticeCore lattice reports
`Periodic()` iff *every* axis wraps, and `Aperiodic()` otherwise.
[`is_bipartite`](@ref) accounts for odd cycles introduced by PBC
with odd-length axes, so the `(-1)^{x+y}` Néel state only reports
as a valid ground state on lattices where it actually is one.

## Further reading

- K. Binder, D. W. Heermann, *Monte Carlo Simulation in Statistical
  Physics*, Ch. 2 (finite-size scaling and BC choice).
- U. Schollwöck, *The density-matrix renormalisation group in the
  age of matrix product states*, Rev. Mod. Phys. **77** (2011).
- A. Gendiar, R. Krcmar, T. Nishino, *Prog. Theor. Phys.* **122**
  (2009) — the original SSD paper.
- N. Byers, C. N. Yang, *Phys. Rev. Lett.* **7** 46 (1961) —
  twisted BC and flux quantisation.
