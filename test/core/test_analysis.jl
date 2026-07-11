using LatticeCore
using LinearAlgebra
using SparseArrays
using Test

# A 1D open chain: its NN graph is a path, so uniform-hopping tight binding
# has the exact particle-in-a-box spectrum -2t·cos(kπ/(N+1)).
chain(N) = LineLattice(N, OpenAxis())

@testset "analysis layer (tight-binding / localization / dynamic SF / regions)" begin
    @testset "tight_binding_hamiltonian reuses adjacency_matrix" begin
        lat = chain(12)
        N = num_sites(lat)
        nb = length(collect(bonds(lat)))
        @test nb == N - 1
        t = 1.3
        H = tight_binding_hamiltonian(lat; t=t)
        # exactly -t on adjacency entries, symmetric, zero diagonal
        @test H ≈ -t .* adjacency_matrix(lat)
        @test issymmetric(Matrix(H))
        @test all(iszero, diag(Matrix(H)))
        @test nnz(H) == 2nb
        @test tr(Matrix(H)^2) ≈ 2nb * t^2 rtol = 1e-12
    end

    @testset "path-graph closed-form spectrum" begin
        for N in (13, 21)
            lat = chain(N)
            t = 0.9
            E = spectrum(lat; t=t)
            exact = sort([-2t * cos(k * π / (N + 1)) for k in 1:N])
            @test E ≈ exact rtol = 1e-10
            @test E ≈ -reverse(E) rtol = 1e-10          # bipartite symmetry
        end
    end

    @testset "on-site shift + per-bond hoppings" begin
        lat = chain(10)
        N = num_sites(lat)
        nb = N - 1
        @test spectrum(lat; onsite=0.7) ≈ spectrum(lat) .+ 0.7 rtol = 1e-10
        @test Matrix(tight_binding_hamiltonian(lat, fill(2.0, nb))) ==
            Matrix(tight_binding_hamiltonian(lat; t=2.0))
        @test_throws DimensionMismatch tight_binding_hamiltonian(lat, fill(1.0, nb + 1))
    end

    @testset "inverse participation ratio" begin
        n = 40
        @test inverse_participation_ratio(fill(1 / sqrt(n), n)) ≈ 1 / n rtol = 1e-12
        e1 = zeros(n);
        e1[3] = 1.0
        @test inverse_participation_ratio(e1) == 1.0
        s = zeros(n);
        s[1:5] .= 1 / sqrt(5)
        @test inverse_participation_ratio(s) ≈ 1 / 5 rtol = 1e-12
        @test_throws ArgumentError inverse_participation_ratio(zeros(n))
    end

    @testset "density of states integrates to N" begin
        lat = chain(21)
        N = num_sites(lat)
        centers, dos = density_of_states(lat; nbins=40)
        @test sum(dos) * (centers[2] - centers[1]) ≈ N rtol = 1e-10
        c2, d2 = density_of_states(lat; nbins=200, broadening=0.05)
        @test sum(d2) * (c2[2] - c2[1]) ≈ N rtol = 5e-2
        @test_throws ArgumentError density_of_states(lat; nbins=0)
    end

    @testset "IPR scaling exponent (extended vs localized)" begin
        Ns = [20, 40, 80, 160, 320]
        @test ipr_scaling_exponent(Ns, [3.0 * n^(-0.7) for n in Ns]) ≈ 0.7 rtol = 1e-10
        @test_throws ArgumentError ipr_scaling_exponent([10], [0.1])
        # extended chain: mean IPR ~ 1/N  ⇒  τ ≈ 1
        chains = [chain(n) for n in (30, 60, 120, 240)]
        s, m = ipr_scaling(chains)
        @test ipr_scaling_exponent(s, m) > 0.9
    end

    @testset "dynamic structure factor A(k,ω): sum rule + first moment" begin
        lat = chain(21)
        N = num_sites(lat)
        t = 1.0
        kpts = [[0.0], [0.5], [1.0], [Float64(π)]]
        omegas, A = dynamic_structure_factor(lat, kpts; t=t, nomega=600, broadening=0.03)
        dω = omegas[2] - omegas[1]
        @test all(A .≥ 0)
        bs = collect(bonds(lat))
        for (j, k) in enumerate(kpts)
            @test sum(@view A[j, :]) * dω ≈ 1.0 rtol = 1e-3     # ∫A dω = ⟨k|k⟩ = 1
            m1 = sum(omegas[m] * A[j, m] for m in eachindex(omegas)) * dω
            exact =
                -2t / N *
                sum(cos(sum(k[d] * b.vector[d] for d in eachindex(k))) for b in bs)
            @test m1 ≈ exact atol = 2e-3                        # ∫ω A dω = ⟨k|H|k⟩
        end
    end

    @testset "2D SimpleSquareLattice tight-binding invariants" begin
        lat = SimpleSquareLattice(5, 5)
        N = num_sites(lat)
        nb = length(collect(bonds(lat)))
        t = 1.0
        E = eigvals(Symmetric(Matrix(tight_binding_hamiltonian(lat; t=t))))
        @test sum(E) ≈ 0 atol = 1e-9
        @test sum(abs2, E) ≈ 2nb * t^2 rtol = 1e-10
        maxdeg = maximum(length(neighbors(lat, i)) for i in 1:N)
        @test maximum(abs, E) ≤ t * maxdeg + 1e-9
    end

    @testset "edge / bulk regions" begin
        lat = chain(10)                       # open ends are the edge
        @test edge_sites(lat) == [1, 10]
        @test edge_sites(lat; depth=2) == [1, 2, 9, 10]
        @test bulk_sites(lat) == collect(2:9)
        eb = edge_bonds(lat)
        @test length(eb) == 2                 # bonds (1-2) and (9-10)
        @test all(b.i ∈ (1, 10) || b.j ∈ (1, 10) for b in eb)
        @test_throws ArgumentError edge_sites(lat; depth=0)
        # periodic chain has no under-coordinated sites
        per = LineLattice(10, PeriodicAxis())
        @test isempty(edge_sites(per))
        @test bulk_sites(per) == collect(1:10)
        @test isempty(edge_bonds(per))
    end
end
