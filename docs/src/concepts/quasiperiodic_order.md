# Quasiperiodic order

This column is still being written. The short version:

## Aperiodic but ordered

Most of condensed matter physics takes periodicity for granted:
a Bravais lattice, translated forever. But there is a family of
structures — Penrose tilings, Ammann–Beenker tilings, Fibonacci
chains, icosahedral quasicrystals — that are *not* periodic in the
usual sense yet are still highly ordered. They have long-range
bond orientational order, crystallographically-forbidden point
symmetries (5-fold for Penrose, 8-fold for Ammann–Beenker, 10-fold
for decagonal quasicrystals), and sharp Bragg peaks in their
diffraction patterns.

## Cut and project

The cleanest way to generate these structures is **cut-and-project**:

1. Start with a periodic lattice in a higher dimension
   ``\mathbb{Z}^{D_{\text{hyper}}}``. For Penrose, ``D_{\text{hyper}} = 5``;
   for Ammann–Beenker, ``4``; for the Fibonacci chain, ``2``.
2. Choose a **physical** subspace of dimension ``D_{\text{phys}}``
   (the plane or line you actually want to look at) and a
   **perpendicular** subspace of dimension
   ``D_{\text{hyper}} - D_{\text{phys}}`` (the "internal" space).
3. Pick an **acceptance window** in the perpendicular space: a
   bounded shape that acts as a stencil.
4. A higher-dimensional lattice point lands in the physical
   structure *if and only if* its perpendicular projection sits
   inside the acceptance window.

The resulting set of physical-space points is the quasicrystal. It
inherits the higher-dimensional translation symmetry but projects
it irrationally, which is why the structure is dense in Fourier
space yet discrete in physical space.

LatticeCore captures the scaffolding with the types
[`AcceptanceWindow`](@ref), [`HyperReciprocalLattice`](@ref), and
[`BraggPeakSet`](@ref); the concrete algorithms that build them
live in `QuasiCrystal.jl`.

## Fourier module and Bragg peaks

Because a quasicrystal inherits its translation symmetry from a
higher-dimensional periodic lattice, its Fourier transform is a
*finite set of Bragg peaks projected through the same matrix*
``\pi_{\parallel}``. Unlike a Bravais reciprocal lattice, the set
of peak locations is **dense** in physical reciprocal space — but
only finitely many peaks are bright, and their intensities are
given by the Fourier transform of the acceptance window evaluated
in the perpendicular space.

A [`BraggPeakSet`](@ref) is what you get after fixing a cutoff:
the finite list of peaks whose intensities exceed some threshold,
stored together with their original higher-dimensional indices so
you can walk back to the hyper lattice.

Crucially, [`BraggPeakSet`](@ref) subtypes
[`AbstractMomentumLattice`](@ref), so the *same*
[`structure_factor`](@ref) code paths that work on a
[`PeriodicMomentumLattice`](@ref) work on a quasicrystal.

## TODO

- Worked example: Fibonacci chain at `depth = 8`, a picture of its
  acceptance window, and the first handful of Bragg peaks.
- Diagram of cut-and-project in 2D projecting down to 1D.
- The role of golden ratio φ in the Fibonacci case and how it shows
  up as the eigenvalue of the substitution matrix.
- When *pseudo-Brillouin zones* are a useful concept for
  quasicrystals and when they break down.

## Further reading

- M. Senechal, *Quasicrystals and Geometry* (Cambridge, 1995).
- M. Baake and U. Grimm, *Aperiodic Order*, Vol. 1 (Cambridge, 2013).
- Z. M. Stadnik (ed.), *Physical Properties of Quasicrystals*
  (Springer, 1999).
