# Reciprocal lattice and Brillouin zone

The real-space picture — points, bonds, sublattices — is half the
story. The other half lives in **momentum space**: the reciprocal
lattice, the Brillouin zone, and the structure factor.
Everything LatticeCore exposes under `reciprocal_lattice`,
`monkhorst_pack`, `BraggPeakSet`, and `structure_factor` ultimately
comes from the content of this column.

## The reciprocal lattice

Given a Bravais lattice with real-space basis
``\mathbf{a}_1, \dots, \mathbf{a}_d``, the **reciprocal basis**
``\mathbf{b}_1, \dots, \mathbf{b}_d`` is defined by the orthogonality
relation

```math
\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi \, \delta_{ij}.
```

Writing ``A`` for the real-space basis matrix and ``B`` for the
reciprocal one,

```math
B = 2\pi \, (A^{-\top}).
```

The **reciprocal lattice** is then the set

```math
\mathbf{G} = m_1 \mathbf{b}_1 + \dots + m_d \mathbf{b}_d,
\qquad m_i \in \mathbb{Z}.
```

Reciprocal lattice vectors are distinguished by a single property:
for every real-space lattice point ``\mathbf{R}``,

```math
e^{-i \mathbf{G} \cdot \mathbf{R}} = 1.
```

In other words, ``\mathbf{G}`` and ``\mathbf{G} + \mathbf{b}_i``
produce *identical* plane waves sampled at the crystal points. Up to
reciprocal lattice translations, only ``\mathbf{k}`` values in a
single primitive cell of reciprocal space are physically distinct.
That primitive cell is the **Brillouin zone**.

In LatticeCore, `reciprocal_basis(ml)` returns the ``B`` matrix for
a [`PeriodicMomentumLattice`](@ref), and the orthogonality relation
is one of the package's unit tests: `a_i · b_j ≈ 2π δ_ij` is
checked directly against the reference square lattice.

## The Brillouin zone

The first Brillouin zone is the Wigner–Seitz cell of the reciprocal
lattice: the set of k-points that are closer to the origin ``\Gamma``
than to any other reciprocal lattice vector. For the square lattice
it is the square ``(-\pi, \pi] \times (-\pi, \pi]``; for the
triangular lattice it is a hexagon; for the honeycomb, a hexagon
rotated by 30°.

The Brillouin zone matters because **every physical quantity that
depends on momentum is periodic in reciprocal space**. Band
structures, structure factors, form factors — they all repeat with
period ``\mathbf{G}``, so all the independent information fits
inside a single BZ.

LatticeCore's reference implementations do not enumerate
Wigner–Seitz cells (that is topology-specific work for
`Lattice2D.jl`). What they *do* provide is a uniform sampler of the
BZ interior — the [`monkhorst_pack`](@ref) mesh.

## Sampling: Monkhorst–Pack and Γ-centred grids

Monte Carlo simulations on a **finite** ``L_1 \times L_2 \times
\dots \times L_d`` lattice can only resolve a finite set of k-points,
namely those compatible with the sample's periodicity:

```math
\mathbf{k} \in \left\{ \frac{n_i}{L_i} \mathbf{b}_i \right\}_{n_i = 0}^{L_i - 1}.
```

LatticeCore ships two equivalent ways to enumerate these points:

- [`gamma_centered(basis, mesh)`](@ref gamma_centered) uses
  ``n_i / L_i``, including the ``\Gamma`` point exactly.
- [`monkhorst_pack(basis, mesh)`](@ref monkhorst_pack) uses
  ``(n_i + 1/2)/L_i - 1/2``, which centres the mesh on ``\Gamma``
  but avoids it. This is the convention solid-state codes use for
  integrating smooth functions over the BZ.

Both return the same number of k-points, ``\prod_i L_i``, matching
the number of sites in the corresponding real-space sample.

## The structure factor

Given a scalar field ``s_i`` defined on the lattice points, the
**structure factor** is the modulus squared of its spatial Fourier
transform:

```math
S(\mathbf{k}) \;=\; \frac{1}{N}\,
    \left| \sum_{i = 1}^{N} s_i \, e^{-i \mathbf{k} \cdot \mathbf{r}_i} \right|^{2}.
```

It tells you, at a single configuration, how much Fourier weight
sits at each momentum. Three sanity checks encode the physics:

- **Ferromagnet**, ``s_i \equiv 1``. The sum is simply ``N``, so
  ``S(0) = N``, and — crucially — ``S(\mathbf{G}) = N`` for *every*
  reciprocal lattice vector, because every such ``\mathbf{G}`` is
  equivalent to ``\Gamma``. Away from ``\mathbf{G}`` the sum over
  sites vanishes by destructive interference. This is verified in
  LatticeCore's test suite.
- **Néel antiferromagnet** on a bipartite lattice,
  ``s_i = (-1)^{\mathbf{x} + \mathbf{y}}`` on a square. The sum
  peaks at ``\mathbf{k} = (\pi, \pi)``, giving ``S(\pi, \pi) = N``
  and ``S(0) = 0``. Bipartiteness is exactly the condition that
  this momentum lies on the reciprocal lattice.
- **Disorder** gives ``S(\mathbf{k}) \approx 1`` on average — a
  useful null baseline when looking for order.

Those checks are what LatticeCore's `test_structure_factor.jl`
actually runs against [`structure_factor`](@ref).

## Implementation notes

LatticeCore uses the **naive O(N) per k-point** formula above
because it is transparent, correct on any lattice — including
quasicrystals, which have no FFT-friendly mesh — and fast enough for
the reference implementations. Production-grade codes on regular
meshes can specialise `structure_factor` for FFT-based evaluation
without changing the API.

For quasicrystals, ``\mathbf{k}`` is no longer drawn from a discrete
reciprocal lattice but from a **dense Fourier module** projected
from a higher-dimensional periodic lattice. LatticeCore ships the
type skeletons [`HyperReciprocalLattice`](@ref) and
[`BraggPeakSet`](@ref); the enumeration algorithm itself lives in
`QuasiCrystal.jl`. See
[Quasiperiodic order](quasiperiodic_order.md) for the column on this.

## Further reading

- Ashcroft–Mermin, Ch. 5 and Ch. 6 (reciprocal lattice and
  Bragg diffraction).
- Monkhorst & Pack, *Phys. Rev. B* **13** 5188 (1976) — the original
  paper on special k-points.
- L. Landau and E. Lifshitz, *Statistical Physics Part 1*, § 116
  (fluctuation–response and the structure factor).
