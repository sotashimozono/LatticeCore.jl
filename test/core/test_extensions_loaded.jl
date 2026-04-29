using LatticeCore
using StaticArrays
using Test

# Smoke test: every weak-dependency extension we ship should be
# discoverable through `Base.get_extension` once the matching package
# is loaded (done in `test/runtests.jl`).
@testset "extensions are loaded" begin
    @test Base.get_extension(LatticeCore, :LatticeCoreFFTWExt) !== nothing
    @test Base.get_extension(LatticeCore, :LatticeCoreNFFTExt) !== nothing
    @test Base.get_extension(LatticeCore, :LatticeCorePlotsExt) !== nothing
end

@testset "NFFT extension installs HasFourierModule fallback" begin
    # The NUFFT path is a stub for now (issue #28). It must still
    # return correct values (matching naive) and not crash. We use
    # a tiny BraggPeakSet as the HasFourierModule lattice.
    bps = BraggPeakSet{1,2,Float64}(
        [SVector(0.0), SVector(1.0), SVector(-1.0)],
        [1.0, 0.3, 0.3],
        [(0, 0), (1, 0), (-1, 0)],
    )
    state = ones(Float64, num_k_points(bps))
    # The stub emits a `@warn` once; suppress with a stderr redirect.
    out = redirect_stderr(devnull) do
        structure_factor(bps, state, bps)
    end
    slow = LatticeCore._structure_factor_naive(bps, state, bps)
    @test out ≈ slow rtol = 1e-12
end
