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
