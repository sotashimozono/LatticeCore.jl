module LatticeCoreNFFTExt

using LatticeCore
using NFFT
using StaticArrays
using LinearAlgebra

"""
    LatticeCoreNFFTExt

Optional `LatticeCore` extension that wires up an NUFFT-based fast path
for `structure_factor(lat, state, ml)` when the lattice carries the
`HasFourierModule` trait — i.e. cut-and-project quasicrystals — by
trait-dispatching on `HasFourierModule` and routing to NFFT.jl.

The NUFFT integration itself is **not yet implemented** at the time
this extension landed (tracked in issue #28). The trait override below
falls back to the naive `O(N · M)` loop with a once-per-session warning,
so the extension is safe to load on top of the existing API: behaviour
is unchanged, only the dispatch site is in place to be filled in later.

Loaded automatically once the user issues `using NFFT`.
"""
LatticeCoreNFFTExt

const _WARNED = Ref(false)

# Trait override for HasFourierModule. For now this delegates to the
# naive helper after emitting a one-shot warning. The intent is to swap
# in NFFT.NNDFTPlan / NFFT.plan_nfft based on whether the supplied
# momentum lattice is regular (`PeriodicMomentumLattice`) or a sparse
# `BraggPeakSet`.
function LatticeCore._structure_factor_fast(
    ::LatticeCore.HasFourierModule,
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice,
)
    if !_WARNED[]
        @warn "LatticeCoreNFFTExt: NUFFT structure_factor path is not yet implemented; falling back to naive O(N · M). See issue #28."
        _WARNED[] = true
    end
    return LatticeCore._structure_factor_naive(lat, state, ml)
end

end # module LatticeCoreNFFTExt
