# Momentum space

LatticeCore's momentum-space layer is designed to be usable
*without* thinking about whether your lattice is periodic or
quasiperiodic. A periodic lattice produces a
[`PeriodicMomentumLattice`](@ref); a quasicrystal produces a
[`BraggPeakSet`](@ref). Both subtype
[`AbstractMomentumLattice`](@ref), and both are themselves
`AbstractLattice{D, T}` — so every observable that walks
`num_sites` / `position` on a real-space lattice transparently
walks `num_k_points` / `k_point` on a momentum lattice.

For the physics background, see
[Reciprocal lattice and Brillouin zone](../concepts/reciprocal_and_brillouin.md).

## The trait-dispatched entry point

Any `AbstractLattice` declares its k-space capability through the
[`reciprocal_support`](@ref) trait:

- [`HasReciprocal`](@ref) — Bravais reciprocal lattice available;
  call [`reciprocal_lattice(lat)`](@ref).
- [`HasFourierModule`](@ref) — a cut-and-project quasicrystal,
  call [`fourier_module(lat)`](@ref) (implemented in
  `QuasiCrystal.jl`).
- [`NoReciprocal`](@ref) — no k-space representation.

The unified helper is [`momentum_lattice(lat)`](@ref), which
dispatches on the trait and throws `ArgumentError` for
`NoReciprocal`:

```julia
lat = SimpleSquareLattice(4, 4, PeriodicAxis())
reciprocal_support(lat)      # HasReciprocal()
ml = momentum_lattice(lat)   # PeriodicMomentumLattice{2, Float64}

obc = SimpleSquareLattice(4, 4, OpenAxis())
reciprocal_support(obc)      # NoReciprocal()
momentum_lattice(obc)        # throws ArgumentError
```

For a periodic lattice, `momentum_lattice` returns the same object
as `reciprocal_lattice(lat)`.

## Building a mesh from a basis

If you know the reciprocal basis and you want to sample the BZ
directly, use [`monkhorst_pack`](@ref) or [`gamma_centered`](@ref).

```julia
using StaticArrays
B = SMatrix{2, 2, Float64}(2π, 0.0,
                           0.0, 2π)

mp = monkhorst_pack(B, (4, 4))          # 16 Monkhorst–Pack k-points
mp isa PeriodicMomentumLattice{2, Float64}

γ = gamma_centered(B, (4, 4))           # 16 Γ-inclusive k-points
```

Both return a `PeriodicMomentumLattice{D, T}` that behaves like
any other lattice: `num_sites(mp)` is 16, `position(mp, 1)` is a
k-vector, and `boundary(mp)` is a uniform periodic wrap (k-space
repeats modulo the reciprocal basis).

The difference between the two is where the mesh sits relative
to Γ:

- `gamma_centered((3, 3))` produces fractions `{0, 1/3, 2/3}` on
  each axis — Γ is one of the sampled points.
- `monkhorst_pack((3, 3))` produces fractions
  `{(n + 0.5)/3 - 0.5} = {-1/3, 0, 1/3}` — the mesh is centred on
  Γ but offset by half a step, which is the convention
  solid-state codes use when integrating smooth functions over
  the BZ.

## Reference-lattice reciprocals

`LineLattice` and `SimpleSquareLattice` both implement
[`basis_vectors`](@ref) (the unit-spacing identity) and
[`reciprocal_lattice`](@ref), which constructs
``B = 2\pi A^{-\top}`` and feeds it through `monkhorst_pack`. On a
square lattice that gives you ``B = 2\pi \, \mathrm{I}_2`` and a
mesh whose density matches the real-space sample:

```julia
lat = SimpleSquareLattice(4, 4, PeriodicAxis())
A = basis_vectors(lat)                   # SMatrix{2, 2, Float64}(1, 0, 0, 1)
ml = reciprocal_lattice(lat)
B = reciprocal_basis(ml)                  # 2π * I

# Orthogonality: a_i · b_j = 2π δ_ij
using LinearAlgebra
[dot(A[:, i], B[:, j]) for i in 1:2, j in 1:2]
# → [2π 0; 0 2π]
```

Non-periodic boundaries (OBC or cylinders) throw
`ArgumentError` from `reciprocal_lattice`, matching the
`reciprocal_support` trait.

## Structure factor

[`structure_factor`](@ref) evaluates

```math
S(\mathbf{k}) = \frac{1}{N}\,
    \left| \sum_i s_i \, e^{-i \mathbf{k} \cdot \mathbf{r}_i} \right|^2
```

for a scalar configuration `state`. The implementation is naive
O(N) per k-point; FFT-based specialisations can be added later
without changing the API.

Three canonical checks:

### Ferromagnet: S(0) = N

```julia
lat = SimpleSquareLattice(4, 4, PeriodicAxis())
state = ones(Int8, num_sites(lat))
structure_factor(lat, state, SVector(0.0, 0.0))   # 16.0
```

### Ferromagnet vanishes at BZ-interior k

```julia
# k = (π/2, 0) is inside the BZ, so Σ_x exp(-i(π/2)x) = 0 for x=1..4
structure_factor(lat, state, SVector(π / 2, 0.0))  # ≈ 0
```

Note that `S` does **not** vanish at `(2π, 0)` — that is a
reciprocal-lattice vector, hence Γ-equivalent, and the
ferromagnetic signal concentrates there too.

### Néel antiferromagnet: S(π, π) = N

```julia
lat = SimpleSquareLattice(4, 4, PeriodicAxis())
N = num_sites(lat)
neel = Int8[
    (-1)^(Int(position(lat, i)[1]) + Int(position(lat, i)[2]))
    for i in 1:N
]

structure_factor(lat, neel, SVector(π, π))    # 16.0
structure_factor(lat, neel, SVector(0.0, 0.0)) # ≈ 0
```

Bipartiteness is exactly the condition that this momentum lies
on the reciprocal lattice of the sample; the 4 × 4 square is
bipartite (both axes even), so the signal is clean.

## Sweeping a whole momentum lattice

You can pass a momentum lattice in place of a single k-vector and
LatticeCore will evaluate the structure factor at every k-point
in the mesh:

```julia
ml = reciprocal_lattice(lat)
Sks = structure_factor(lat, neel, ml)        # Vector{Float64}, length num_sites(ml)
argmax(Sks)                                   # the k-index where the peak lives
```

The result is a plain `Vector{Float64}` of the same length as
`num_k_points(ml)`, indexed the same way as `k_point(ml, i)`.

## Quasicrystal Fourier modules (preview)

Everything above generalises: a [`BraggPeakSet`](@ref) subtypes
[`AbstractMomentumLattice`](@ref), exposes
`num_k_points(bps) = length(bps.peaks)`, and answers
`position(bps, i) = bps.peaks[i]`. That means the *same*
`structure_factor(lat, state, ml)` code path works for a
quasicrystal — there is no "quasi" branch in the observer. You
iterate over the Bragg peaks and each one returns a plain
`Float64`.

The generation algorithm — choosing an acceptance window, walking
the higher-dimensional reciprocal lattice, computing the window
Fourier transform — lives in `QuasiCrystal.jl`, not here.
LatticeCore ships the type skeletons
[`HyperReciprocalLattice`](@ref) and [`AcceptanceWindow`](@ref)
so downstream packages have a stable interface to target.

## See also

- [Concepts: Reciprocal lattice and Brillouin zone](../concepts/reciprocal_and_brillouin.md)
- [Concepts: Quasiperiodic order](../concepts/quasiperiodic_order.md)
- [Reference: Momentum space](../reference/momentum.md)
