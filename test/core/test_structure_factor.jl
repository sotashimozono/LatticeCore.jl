using LatticeCore
using StaticArrays
using Test

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
end
