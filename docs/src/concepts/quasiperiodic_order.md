# Quasiperiodic order

Most of condensed-matter physics takes **periodicity** for granted.
Crystals are invariant under a Bravais lattice of translations, and
almost every technique — Bloch theorem, band structure, reciprocal
lattices — is a consequence. But in 1982 Shechtman's diffraction
measurements of a rapidly cooled Al–Mn alloy showed a pattern with
sharp Bragg peaks *and* 10-fold rotational symmetry, which is
forbidden for any 2D or 3D crystallographic space group. The
discovery of **quasicrystals** (Nobel Prize in Chemistry, 2011)
forced physics to accept a structural regime that is neither
periodic nor disordered: long-range order without translation
symmetry.

LatticeCore is designed to host these structures on the same
footing as Bravais crystals. This column is the physics story; the
construction and Fourier analysis algorithms live in
`QuasiCrystal.jl` and are exposed here through the type skeletons
[`HyperReciprocalLattice`](@ref), [`AcceptanceWindow`](@ref), and
[`BraggPeakSet`](@ref).

## What "quasiperiodic" means

A point set ``\mathcal{L} \subset \mathbb{R}^d`` is **quasiperiodic**
if its autocorrelation function — equivalently, its diffraction
pattern — is a **pure point measure** supported on a *dense*
countable set:

```math
\mathcal{F}[\rho](\mathbf{k}) = \sum_{\mathbf{k}_n \in M^\ast}
    I(\mathbf{k}_n)\, \delta(\mathbf{k} - \mathbf{k}_n).
```

Periodic crystals are a special case where the support
``M^\ast`` is a discrete Bravais reciprocal lattice. For a
quasicrystal, ``M^\ast`` is still countable but **dense** in
``\mathbb{R}^d``: you can find a Bragg peak arbitrarily close to
any ``\mathbf{k}``, yet only countably many peaks carry non-zero
intensity. The resulting diffraction pattern is sharp but
inexhaustible.

The intensities ``I(\mathbf{k}_n)`` are not uniform. They fall off
with a form factor, leaving a finite number of *bright* peaks at
any intensity threshold — which is exactly the finite set
LatticeCore stores in a [`BraggPeakSet`](@ref).

## Cut and project

The cleanest way to generate a quasicrystal is **cut-and-project**.
The recipe:

1. Pick a host **periodic** lattice in a higher dimension
   ``D_{\text{hyper}}``. For the Fibonacci chain, ``D_{\text{hyper}} = 2``;
   for the Penrose tiling, ``D_{\text{hyper}} = 5``; for the
   Ammann–Beenker tiling, ``D_{\text{hyper}} = 4``.
2. Choose a **physical subspace** ``E_\parallel`` of dimension
   ``D_{\text{phys}}`` inside ``\mathbb{R}^{D_{\text{hyper}}}``.
   The embedding must be irrational — that is, ``E_\parallel`` is
   not spanned by any rational basis of the hyper lattice — so
   that its parallel projection lands *densely* in ``E_\parallel``
   without ever repeating.
3. Choose a **perpendicular subspace** ``E_\perp`` of dimension
   ``D_{\text{hyper}} - D_{\text{phys}}``, orthogonal to
   ``E_\parallel``.
4. Choose an **acceptance window** ``W \subset E_\perp`` — a
   bounded region that acts as a stencil.
5. Project every hyper lattice point ``\mathbf{n} \in
   \mathbb{Z}^{D_{\text{hyper}}}`` whose **perpendicular**
   projection sits inside ``W`` onto the physical subspace. The
   resulting set of physical-space points

```math
   \mathcal{L} = \{\, \pi_\parallel(\mathbf{n}) :
                      \mathbf{n} \in \mathbb{Z}^{D_{\text{hyper}}},
                      \pi_\perp(\mathbf{n}) \in W \,\}
```

is the quasicrystal.

The acceptance window is where all the structure hides. A
rectangular window in 2D gives you a Fibonacci chain. Five
pentagonal windows arranged with 5-fold symmetry (one per
perpendicular-space slab of the Penrose lifting) give you the
Penrose P3 tiling. Octagonal windows give Ammann–Beenker.

LatticeCore captures the three ingredients as type skeletons:

```julia
abstract type AcceptanceWindow end

struct HyperReciprocalLattice{DPhys, DHyper, DPerp, T, W <: AcceptanceWindow}
    hyper_basis::SMatrix{DHyper, DHyper, T}
    parallel_proj::SMatrix{DPhys, DHyper, T}
    perp_proj::SMatrix{DPerp, DHyper, T}
    window::W
end
```

Concrete window types (`PolygonalWindow`, `IntervalWindow`, ...)
and the construction of projection matrices live in
`QuasiCrystal.jl`. The shape lives here because it is the contract
every downstream consumer agrees on.

## Fibonacci in full

The Fibonacci chain is the simplest non-trivial example, small
enough to do by hand:

- Hyper lattice: ``\mathbb{Z}^2``.
- Physical subspace: the line with slope ``1/\varphi`` where
  ``\varphi = (1 + \sqrt{5})/2`` is the golden ratio.
- Perpendicular subspace: the orthogonal line with slope
  ``-\varphi``.
- Acceptance window: an interval of length ``\sqrt{1 + \varphi^2}``
  on the perpendicular line.

Project every integer point ``(m, n) \in \mathbb{Z}^2`` whose
perpendicular projection lies in the window, and you get a
sequence of "long" (``L``) and "short" (``S``) intervals on the
physical line. Asymptotically the ratio ``|L|/|S| = \varphi``, and
the substitution rule ``L \to LS``, ``S \to L`` iterated from a
single ``L`` generates the same sequence.

The substitution matrix

```math
M = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}
```

has eigenvalues ``\varphi`` and ``-1/\varphi``. The leading
eigenvalue ``\varphi`` is the **inflation factor**: every
substitution step stretches the sequence by ``\varphi``. This
simple fact — an algebraic integer lurking in a hyperbolic map —
is what makes Fibonacci quasicrystals so tractable analytically.

## The Fourier module and Bragg peaks

Because a quasicrystal inherits its translation symmetry from the
host periodic lattice, its Fourier transform is supported on the
**projected reciprocal lattice** of the host:

```math
M^\ast = \pi_\parallel (\mathbb{Z}^{D_{\text{hyper}}})^\ast.
```

This projected set is countable (it is the image of a countable
set) and, crucially, **dense** in ``\mathbb{R}^{D_{\text{phys}}}``
— because the projection is irrational. The set ``M^\ast`` is the
**Fourier module** of the quasicrystal. Unlike a Bravais
reciprocal lattice, ``M^\ast`` is not discrete.

The intensity carried by each projected peak is the Fourier
transform of the acceptance window evaluated in the perpendicular
space:

```math
I(\mathbf{k}_n) \propto \bigl|\, \hat{W}(\pi_\perp(\mathbf{g}_n))\, \bigr|^2,
```

where ``\mathbf{g}_n \in (\mathbb{Z}^{D_{\text{hyper}}})^\ast`` is
the higher-dimensional reciprocal lattice vector that
``\mathbf{k}_n = \pi_\parallel(\mathbf{g}_n)`` projects from. For
an interval window of width ``a``,
``\hat{W}(q) \propto \mathrm{sinc}(qa/2)``; for a polygonal
window, the window Fourier transform is a sum of sinc factors
coming from the polygon's vertices.

The practical consequence is that **only a finite number of
peaks exceed any given intensity threshold**. Pick a cutoff, list
those peaks, and you have a [`BraggPeakSet`](@ref):

```julia
struct BraggPeakSet{DPhys, DHyper, T} <: AbstractMomentumLattice{DPhys, T}
    peaks::Vector{SVector{DPhys, T}}
    intensities::Vector{T}
    hyper_indices::Vector{NTuple{DHyper, Int}}
end
```

This is a subtype of [`AbstractMomentumLattice`](@ref), so the
*same* [`structure_factor`](@ref) code paths that work on a
[`PeriodicMomentumLattice`](@ref) work on a quasicrystal. In the
words of the 05 momentum-space chapter, a Bragg peak set is the
quasicrystal's answer to a Monkhorst–Pack mesh.

## Why this matters for LatticeCore

You could, in principle, stick a quasicrystal generator on its own
hierarchy and never share interfaces with periodic lattices. That
is what most of the existing quasicrystal literature does. The
argument for the LatticeCore approach — one `AbstractLattice`
hierarchy, one `AbstractMomentumLattice` hierarchy, one structure
factor — is that **you want to run the same Monte Carlo on both**.
A diluted antiferromagnet on a Penrose tiling wants the same
observer stack as the same model on a square lattice: the trait
machinery picks the right reciprocal-space object through
[`momentum_lattice`](@ref), the structure factor iterates over it
with exactly the same code, and the physics stays comparable.

The cut-and-project recipe is how the hierarchy is populated; the
interface is what makes the comparison cheap.

## Further reading

- D. Shechtman *et al.*, *Phys. Rev. Lett.* **53** 1951 (1984) —
  the paper that identified icosahedral quasicrystals.
- M. Senechal, *Quasicrystals and Geometry*, Cambridge (1995).
- M. Baake, U. Grimm, *Aperiodic Order*, Vol. 1, Cambridge (2013)
  — the definitive reference on cut-and-project constructions and
  their Fourier analysis.
- N. G. de Bruijn, *Algebraic theory of Penrose's non-periodic
  tilings*, Kon. Neder. Akad. Wetensch. Proc. Ser. A **84** (1981)
  — the original 5D lift of the Penrose tiling.
- Z. M. Stadnik (ed.), *Physical Properties of Quasicrystals*,
  Springer (1999) — the experimental side.
