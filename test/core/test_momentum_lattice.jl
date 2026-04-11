using LatticeCore
using StaticArrays
using LinearAlgebra
using Test

@testset "AbstractMomentumLattice" begin
    @testset "type hierarchy" begin
        @test AbstractMomentumLattice <: AbstractLattice
        @test PeriodicMomentumLattice <: AbstractMomentumLattice
        @test BraggPeakSet <: AbstractMomentumLattice
    end

    @testset "monkhorst_pack 1D" begin
        basis = SMatrix{1,1,Float64}(2π)          # reciprocal of unit-spacing chain
        ml = monkhorst_pack(basis, (4,))
        @test ml isa PeriodicMomentumLattice{1,Float64}
        @test num_k_points(ml) == 4
        @test num_sites(ml) == 4                  # AbstractLattice delegation
        @test dimension(ml) == 1
        @test reciprocal_basis(ml) == basis

        # MP mesh is half-shifted and Γ-centred: for N=4 the fractions
        # are (n + 0.5)/N - 0.5 = -3/8, -1/8, 1/8, 3/8
        expected_fracs = [-3 / 8, -1 / 8, 1 / 8, 3 / 8]
        for i in 1:4
            @test k_point(ml, i)[1] ≈ 2π * expected_fracs[i]
        end

        @test position(ml, 2) == k_point(ml, 2)   # AbstractLattice delegation
        @test size_trait(ml) isa FiniteSize{1}
        @test is_finite(ml) == true
        @test neighbors(ml, 1) == Int[]           # graph interface is vacuous
    end

    @testset "gamma_centered 2D" begin
        basis = SMatrix{2,2,Float64}(2π, 0.0, 0.0, 2π)
        ml = gamma_centered(basis, (3, 3))
        @test num_k_points(ml) == 9
        # First k-point should be Γ = (0, 0)
        @test k_point(ml, 1) == SVector(0.0, 0.0)
    end

    @testset "BraggPeakSet skeleton (QuasiCrystal will fill it)" begin
        peaks = [SVector(1.0, 2.0), SVector(3.0, 4.0)]
        intensities = [0.5, 0.2]
        hyper_indices = [(1, 0, 0, 1, 0), (0, 1, 1, 0, 0)]   # 5D indices
        bps = BraggPeakSet{2,5,Float64}(peaks, intensities, hyper_indices)
        @test bps isa AbstractMomentumLattice{2,Float64}
        @test num_k_points(bps) == 2
        @test k_point(bps, 1) == SVector(1.0, 2.0)
        @test reciprocal_support(bps) isa HasFourierModule
    end
end

@testset "Trait-dispatched momentum_lattice entry point" begin
    # LineLattice under PBC: HasReciprocal, returns PeriodicMomentumLattice
    lat_pbc = LineLattice(4, PeriodicAxis())
    @test reciprocal_support(lat_pbc) isa HasReciprocal
    ml = momentum_lattice(lat_pbc)
    @test ml isa PeriodicMomentumLattice{1,Float64}
    @test num_k_points(ml) == 4

    # LineLattice under OBC: NoReciprocal, momentum_lattice throws
    lat_obc = LineLattice(4, OpenAxis())
    @test reciprocal_support(lat_obc) isa NoReciprocal
    @test_throws ArgumentError momentum_lattice(lat_obc)
end

@testset "LineLattice reciprocal_lattice" begin
    lat = LineLattice(6, PeriodicAxis())
    @test basis_vectors(lat) == SMatrix{1,1,Float64}(1.0)

    ml = reciprocal_lattice(lat)
    @test ml isa PeriodicMomentumLattice{1,Float64}
    @test reciprocal_basis(ml)[1, 1] ≈ 2π
    @test num_k_points(ml) == 6

    # OBC → should error
    @test_throws ArgumentError reciprocal_lattice(LineLattice(6, OpenAxis()))
end

@testset "SimpleSquareLattice reciprocal_lattice" begin
    lat = SimpleSquareLattice(3, 3, PeriodicAxis())
    A = basis_vectors(lat)
    @test A == SMatrix{2,2,Float64}(1.0, 0.0, 0.0, 1.0)

    ml = reciprocal_lattice(lat)
    @test ml isa PeriodicMomentumLattice{2,Float64}
    B = reciprocal_basis(ml)
    # Orthogonality: a_i · b_j = 2π δ_ij
    for i in 1:2, j in 1:2
        expected = i == j ? 2π : 0.0
        @test dot(A[:, i], B[:, j]) ≈ expected
    end

    @test num_k_points(ml) == 9

    # Cylinder / OBC → error
    cyl = SimpleSquareLattice(3, 3, LatticeBoundary((PeriodicAxis(), OpenAxis())))
    @test_throws ArgumentError reciprocal_lattice(cyl)
end
