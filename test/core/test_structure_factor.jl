using LatticeCore
using StaticArrays
using Test

# Dummy wrapper used to verify that `LatticeCore._has_known_grid` /
# `LatticeCore._reshape_state` can be overridden by a *downstream*
# Bravais lattice without depending on FFTW. Defined at top level so
# the type is reachable from method tables (a `@testset` block opens
# a `let`-style local scope, which would prevent that).
struct DummyBravais{L} <: AbstractLattice{2,Float64}
    inner::L
end
LatticeCore.position(d::DummyBravais, i::Int) = position(d.inner, i)
LatticeCore.num_sites(d::DummyBravais) = num_sites(d.inner)
LatticeCore.neighbors(d::DummyBravais, i::Int) = neighbors(d.inner, i)
LatticeCore.boundary(d::DummyBravais) = boundary(d.inner)
LatticeCore.size_trait(d::DummyBravais) = LatticeCore.size_trait(d.inner)
LatticeCore.reciprocal_support(d::DummyBravais) = LatticeCore.reciprocal_support(d.inner)
LatticeCore.site_layout(d::DummyBravais) = LatticeCore.site_layout(d.inner)

@testset "structure_factor" begin
    @testset "ferromagnet at k=0" begin
        # S(k=0) = |Σ s_i|² / N = N² / N = N for s_i ≡ 1
        lat = SimpleSquareLattice(4, 4, PeriodicAxis())
        state = ones(Int8, num_sites(lat))
        S0 = structure_factor(lat, state, SVector(0.0, 0.0))
        @test S0 ≈ Float64(num_sites(lat))
    end

    @testset "ferromagnet vanishes at BZ-interior k" begin
        # A k inside the first BZ but not a reciprocal-lattice vector.
        # At k = (π/2, 0), Σ_{x=1..4} exp(-i(π/2)x) = -i + (-1) + i + 1 = 0,
        # so the row sum vanishes and S(k) = 0.
        lat = SimpleSquareLattice(4, 4, PeriodicAxis())
        state = ones(Int8, num_sites(lat))
        S_mid = structure_factor(lat, state, SVector(π / 2, 0.0))
        @test isapprox(S_mid, 0.0; atol=1e-10)

        # At a full reciprocal-lattice vector b = (2π, 0), the phase is
        # trivial on every integer site, so we are still sitting at Γ:
        # S(b) = N, matching S(0). Sanity-check this to lock the
        # convention down.
        S_G = structure_factor(lat, state, SVector(2π, 0.0))
        @test S_G ≈ Float64(num_sites(lat))
    end

    @testset "Neel antiferromagnet peaks at (π, π)" begin
        # 4x4 PBC square is bipartite; build the Neel state s = (-1)^(x+y).
        lat = SimpleSquareLattice(4, 4, PeriodicAxis())
        N = num_sites(lat)
        state = Vector{Int8}(undef, N)
        for i in 1:N
            p = position(lat, i)
            state[i] = Int8((-1)^(Int(p[1]) + Int(p[2])))
        end
        # S(π, π) should equal N for the Neel state.
        S_afm = structure_factor(lat, state, SVector(Float64(π), Float64(π)))
        @test S_afm ≈ Float64(N)

        # Off-peak should be small at (0, 0)
        S0 = structure_factor(lat, state, SVector(0.0, 0.0))
        @test isapprox(S0, 0.0; atol=1e-10)
    end

    @testset "structure_factor over a whole momentum lattice" begin
        lat = SimpleSquareLattice(4, 4, PeriodicAxis())
        state = ones(Int8, num_sites(lat))
        ml = reciprocal_lattice(lat)
        Sks = structure_factor(lat, state, ml)
        @test Sks isa Vector{Float64}
        @test length(Sks) == num_k_points(ml)
        # Ferromagnet: all k-space mass is at k = 0. The MP mesh on a
        # 4x4 square does not contain (0, 0), but the total weight
        # should still concentrate on small-k points. A weaker check:
        # the total sum of |S(k)| * (peak indicator) is finite.
        @test all(S >= -1e-10 for S in Sks)
    end

    @testset "1D chain" begin
        lat = LineLattice(6, PeriodicAxis())
        state = ones(Int8, num_sites(lat))
        # k = 0: all phases are 1 → |N|² / N = N
        @test structure_factor(lat, state, SVector(0.0)) ≈ Float64(num_sites(lat))
    end

    # ---- FFT fast path (Bravais, HasReciprocal) ---------------------
    #
    # `LatticeCoreFFTWExt` is loaded eagerly in `test/runtests.jl`, so
    # `structure_factor(lat, state, ml)` should hit the FFT path on
    # `SimpleSquareLattice` / `LineLattice`. The result must match the
    # naive double loop (modulo FP noise).

    @testset "FFT path matches naive on SimpleSquareLattice (MP mesh)" begin
        for (Lx, Ly) in ((4, 4), (6, 4), (8, 8))
            lat = SimpleSquareLattice(Lx, Ly, PeriodicAxis())
            ml = reciprocal_lattice(lat)
            for trial in 1:3
                state = randn(num_sites(lat))
                fast = structure_factor(lat, state, ml)
                slow = LatticeCore._structure_factor_naive(lat, state, ml)
                @test fast ≈ slow rtol = 1e-9
            end
        end
    end

    @testset "FFT path matches naive on LineLattice (MP mesh)" begin
        for N in (4, 6, 16, 32)
            lat = LineLattice(N, PeriodicAxis())
            for trial in 1:3
                state = randn(num_sites(lat))
                ml = reciprocal_lattice(lat)
                fast = structure_factor(lat, state, ml)
                slow = LatticeCore._structure_factor_naive(lat, state, ml)
                @test fast ≈ slow rtol = 1e-9
            end
        end
    end

    @testset "FFT path matches naive on a Γ-centred mesh" begin
        # `gamma_centered` uses zero offset, exercising the
        # no-prephase branch in the FFT extension.
        lat = SimpleSquareLattice(4, 4, PeriodicAxis())
        A = LatticeCore.basis_vectors(lat)
        B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
        ml = LatticeCore.gamma_centered(B, (lat.Lx, lat.Ly))
        state = randn(num_sites(lat))
        fast = structure_factor(lat, state, ml)
        slow = LatticeCore._structure_factor_naive(lat, state, ml)
        @test fast ≈ slow rtol = 1e-9
    end

    @testset "FFT path falls back to naive when mesh ≠ lattice dims" begin
        lat = SimpleSquareLattice(4, 4, PeriodicAxis())
        # A coarser MP mesh — extension should detect the mismatch and
        # fall back to the naive helper.
        A = LatticeCore.basis_vectors(lat)
        B = SMatrix{2,2,Float64}(2π * inv(transpose(A)))
        ml_coarse = LatticeCore.monkhorst_pack(B, (2, 2))
        state = randn(num_sites(lat))
        out = structure_factor(lat, state, ml_coarse)
        slow = LatticeCore._structure_factor_naive(lat, state, ml_coarse)
        @test out ≈ slow rtol = 1e-12
        @test length(out) == 4
    end

    # ---- Downstream opt-in hook -------------------------------------
    #
    # The grid-layout hooks `LatticeCore._has_known_grid` /
    # `LatticeCore._reshape_state` live on `LatticeCore` itself so a
    # downstream package can register a custom Bravais lattice on the
    # FFT fast path without depending on FFTW. Verify the dispatch
    # path works end-to-end by wiring up a tiny clone of
    # `SimpleSquareLattice`.

    @testset "downstream lattices can opt into the FFT fast path" begin
        # Reference lattices ship with the opt-in baked in.
        @test LatticeCore._has_known_grid(SimpleSquareLattice(2, 2, PeriodicAxis())) == true
        @test LatticeCore._has_known_grid(LineLattice(4, PeriodicAxis())) == true

        inner = SimpleSquareLattice(4, 4, PeriodicAxis())
        dlat = DummyBravais(inner)

        # Default: dummy is not on the opt-in list. Confirm before
        # we exercise the override.
        @test LatticeCore._has_known_grid(dlat) == false
        @test_throws ErrorException LatticeCore._reshape_state(dlat, ones(16), (4, 4))

        # Without opt-in the FFT extension's `_fft_eligible` returns
        # `false`, so `structure_factor(dlat, state, ml)` runs the
        # naive helper. The numbers must agree with the reference.
        ml = reciprocal_lattice(inner)
        state = randn(num_sites(inner))
        baseline = LatticeCore._structure_factor_naive(inner, state, ml)
        out_default = structure_factor(dlat, state, ml)
        @test out_default ≈ baseline rtol = 1e-12

        # Now opt the dummy into the fast path. After the override,
        # `_has_known_grid` returns `true` and `_reshape_state` no
        # longer errors. With FFTW loaded eagerly in `runtests.jl`,
        # `structure_factor(dlat, state, ml)` should still match the
        # naive baseline (just routed via the FFT path instead).
        LatticeCore._has_known_grid(::DummyBravais) = true
        LatticeCore._reshape_state(::DummyBravais, s, dims) = reshape(s, dims)

        @test LatticeCore._has_known_grid(dlat) == true
        @test LatticeCore._reshape_state(dlat, state, (4, 4)) isa AbstractMatrix{Float64}

        out_opted_in = structure_factor(dlat, state, ml)
        @test out_opted_in ≈ baseline rtol = 1e-9
    end
end
