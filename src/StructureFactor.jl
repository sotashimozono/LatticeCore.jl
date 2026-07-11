"""
    structure_factor(lat, state, k::SVector) → Float64

Scalar structure factor at a single k-point:

    S(k) = (1/N) |⟨ Σ_i s_i exp(-i k · r_i) ⟩|²

computed on a single snapshot (no thermal averaging). The default
implementation is naive O(N) per k-point; FFT- and NUFFT-based
specialisations on regular meshes are provided in the optional
`LatticeCoreFFTWExt` / `LatticeCoreNFFTExt` extensions and dispatch
on `reciprocal_support(lat)`.
"""
function structure_factor(
    lat::AbstractLattice, state::AbstractVector, k::SVector{D,T}
) where {D,T}
    N = num_sites(lat)
    acc = zero(ComplexF64)
    for i in 1:N
        acc += state[i] * cis(-dot(k, position(lat, i)))
    end
    return abs2(acc) / N
end

"""
    structure_factor(lat, state, ml::AbstractMomentumLattice) → Vector{Float64}

Evaluate the structure factor at every k-point of `ml`.

Dispatch is driven by `reciprocal_support(lat)`:

- `HasReciprocal()` (Bravais periodic lattices): uses an FFT-based
  fast path when `LatticeCoreFFTWExt` is loaded (`using FFTW`),
  otherwise falls back to the naive `O(N · M)` loop.
- `HasFourierModule()` (cut-and-project quasicrystals): uses a NUFFT
  fast path when `LatticeCoreNFFTExt` is loaded (`using NFFT`),
  otherwise falls back to naive.
- `NoReciprocal()` and any unhandled case: naive.

Extensions hook into the `_structure_factor_fast` indirection below;
the default fallback simply runs the naive loop.
"""
function structure_factor(
    lat::AbstractLattice, state::AbstractVector, ml::AbstractMomentumLattice
)
    return _structure_factor_fast(reciprocal_support(lat), lat, state, ml)
end

# Default fallback: naive double loop. Extensions can specialise on the
# trait type (first argument) to install a faster path.
function _structure_factor_fast(
    ::AbstractReciprocalSupport,
    lat::AbstractLattice,
    state::AbstractVector,
    ml::AbstractMomentumLattice,
)
    return _structure_factor_naive(lat, state, ml)
end

"""
    LatticeCore._structure_factor_naive(lat, state, ml) → Vector{Float64}

Reference O(N · M) loop. Kept as the unconditional fallback used by
extensions when the regular-mesh / NUFFT preconditions don't hold.
"""
function _structure_factor_naive(
    lat::AbstractLattice, state::AbstractVector, ml::AbstractMomentumLattice
)
    out = Vector{Float64}(undef, num_k_points(ml))
    for i in 1:num_k_points(ml)
        out[i] = structure_factor(lat, state, k_point(ml, i))
    end
    return out
end

# ---- FFT fast-path opt-in hooks --------------------------------------
#
# These two helpers describe a Bravais lattice's "natural grid layout":
# whether `state[i]` can be reshaped onto a regular `D`-dimensional
# mesh that lines up with the lattice's site index ordering. They are
# defined here (rather than inside `LatticeCoreFFTWExt`) so that
# downstream packages — `Lattice2D`, `QuasiCrystal`, custom lattices —
# can opt their own lattices into the FFT fast path **without depending
# on FFTW**:
#
# ```julia
# # In a downstream package, after `using LatticeCore`:
# struct MyLattice <: AbstractLattice{2,Float64} ... end
# LatticeCore._has_known_grid(::MyLattice) = true
# LatticeCore._reshape_state(::MyLattice, state, dims) = reshape(state, dims)
# ```
#
# When `LatticeCoreFFTWExt` is loaded (`using FFTW`), it dispatches
# through these names; when it isn't, the hooks simply have no effect.
# The defaults are deliberately conservative: an unknown lattice is not
# eligible, and `_reshape_state` errors if called without an opt-in
# (the FFT extension only reaches it after `_has_known_grid` returns
# `true`, so the error path is purely defensive).

"""
    LatticeCore._has_known_grid(lat::AbstractLattice) → Bool

Whether the lattice has a known, FFT-compatible grid layout that
matches a `PeriodicMomentumLattice` mesh of the same `dims`.

Default: `false`. Concrete lattices opt in by adding a method that
returns `true`, paired with a [`_reshape_state`](@ref
LatticeCore._reshape_state) method that produces the corresponding
`D`-dimensional grid view.

The FFT fast path in `LatticeCoreFFTWExt` only fires when this
returns `true`; lattices without an opt-in stay on the naive helper.
"""
_has_known_grid(::AbstractLattice) = false

"""
    LatticeCore._reshape_state(lat, state, dims) → AbstractArray

Reshape a per-site `state::AbstractVector` into the lattice's
natural `D`-dimensional grid of size `dims`, so that the FFT
extension can take a single `fft` over the result.

Default: throws — concrete lattices that opt into
[`_has_known_grid`](@ref LatticeCore._has_known_grid) must also
provide a method here.
"""
function _reshape_state(lat::AbstractLattice, state, dims)
    return error(
        "LatticeCore._reshape_state not defined for $(typeof(lat)); ",
        "opt in by defining `LatticeCore._reshape_state(::MyLattice, state, dims)` ",
        "alongside `LatticeCore._has_known_grid(::MyLattice) = true`.",
    )
end
