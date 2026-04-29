module LatticeCoreFFTWExt

using LatticeCore
using FFTW
using StaticArrays
using LinearAlgebra

"""
    LatticeCoreFFTWExt

Optional `LatticeCore` extension that installs an FFT-based fast path
for `structure_factor(lat, state, ml)` when the lattice carries the
`HasReciprocal` trait and the supplied momentum lattice is a regular
mesh whose dimensions match the real-space lattice grid.

Loaded automatically once the user issues `using FFTW`. Falls back to
the naive double loop whenever the alignment preconditions don't hold,
so adding `using FFTW` never changes results — only performance.

Currently the FFT path is wired up for the reference lattices shipped
with LatticeCore (`LineLattice`, `SimpleSquareLattice`). Downstream
Bravais lattices stay on the naive path until they opt in.
"""
LatticeCoreFFTWExt

# ---- Trait override --------------------------------------------------

# Specialise on `HasReciprocal`. The default fallback in
# `src/StructureFactor.jl` runs the naive loop; we override here when
# the FFTW extension is active, falling back to the naive helper if the
# concrete lattice isn't on the FFT-eligible list.
function LatticeCore._structure_factor_fast(
    ::LatticeCore.HasReciprocal,
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
)
    if _fft_eligible(lat, ml)
        return _structure_factor_fft(lat, state, ml)
    else
        return LatticeCore._structure_factor_naive(lat, state, ml)
    end
end

# ---- Eligibility -----------------------------------------------------

# FFT path requires:
# - ml is a PeriodicMomentumLattice (regular k-grid)
# - lat is a finite Bravais lattice with grid dims matching ml.mesh
# - lat is one of the lattices for which we know the natural grid layout
#
# The grid-layout opt-in lives on `LatticeCore` itself
# (`_has_known_grid` / `_reshape_state`) so downstream packages such as
# `Lattice2D` and `QuasiCrystal` can register their lattices on the
# fast path without taking an FFTW dependency.
function _fft_eligible(
    lat::LatticeCore.AbstractLattice, ml::LatticeCore.AbstractMomentumLattice
)
    ml isa LatticeCore.PeriodicMomentumLattice || return false
    LatticeCore.is_finite(lat) || return false
    st = LatticeCore.size_trait(lat)
    st isa LatticeCore.FiniteSize || return false
    st.dims == ml.mesh || return false
    return LatticeCore._has_known_grid(lat)
end

# ---- FFT core --------------------------------------------------------

# Compute S(k) at every k-point of `ml` via a single FFT.
#
# For a Bravais lattice with `r = origin + Σ_d a_d * n_d` and a regular
# mesh `k = B * frac` with `frac_d = (idx_d - 1 + offset_d) / N_d` and
# `B^T A = 2π I` (the standard Bravais ↔ reciprocal pairing) we have
#
#     k · r = (B · frac) · origin + 2π Σ_d frac_d * n_d
#           = (B · frac) · origin
#             + 2π Σ_d (offset_d / N_d) * n_d
#             + 2π Σ_d (idx_d - 1) * n_d / N_d.
#
# The first two terms are independent of `n` and `idx` respectively — or
# in the case of the origin term, depend only on `idx`. After taking
# `|·|² / N` to form `S(k)`, the per-`idx` phase drops out, leaving a
# DFT of `state .* prephase`, where
#
#     prephase[n] = exp(-i 2π Σ_d (offset_d / N_d) n_d).
#
# The mesh offset is recovered from `B \ k_point(ml, 1)` (which equals
# the fractional coordinate of the first k-point); `idx = 1` corresponds
# to `n = 0` per axis so that fractional coordinate is exactly
# `offset / N`.
function _structure_factor_fft(
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector,
    ml::LatticeCore.PeriodicMomentumLattice{D,T},
) where {D,T}
    dims = ml.mesh
    N = prod(dims)
    grid = LatticeCore._reshape_state(lat, state, dims)

    B = LatticeCore.reciprocal_basis(ml)
    # Fractional offset of the mesh (per axis, scaled by 1/N_d).
    frac0 = SVector{D,Float64}(B \ LatticeCore.k_point(ml, 1))

    grid_c = _apply_prephase(grid, frac0)

    F = fft(grid_c)
    out = Vector{Float64}(undef, N)
    @inbounds for i in eachindex(F)
        out[i] = abs2(F[i]) / N
    end
    return out
end

# Build the per-site complex array `state[n] * exp(-i 2π Σ_d frac0_d * n_d)`.
# When the offset is exactly zero (Γ-centred mesh) we skip the phase
# multiplication entirely.
function _apply_prephase(grid::AbstractArray{<:Any,D}, frac0::SVector{D,Float64}) where {D}
    if all(iszero, frac0)
        return ComplexF64.(grid)
    end
    out = Array{ComplexF64,D}(undef, size(grid))
    @inbounds for I in CartesianIndices(grid)
        s = 0.0
        for d in 1:D
            s += frac0[d] * (I[d] - 1)
        end
        out[I] = ComplexF64(grid[I]) * cis(-2π * s)
    end
    return out
end

end # module LatticeCoreFFTWExt
