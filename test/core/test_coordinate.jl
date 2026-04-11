using LatticeCore
using StaticArrays
using Test

@testset "AbstractCoordinate" begin
    @testset "type hierarchy" begin
        @test RealSpace <: AbstractCoordinate
        @test LatticeCoord <: AbstractCoordinate
        @test HigherDimCoord <: AbstractCoordinate
    end

    @testset "RealSpace construction" begin
        rs_sv = RealSpace(SVector(1.0, 2.0))
        @test rs_sv isa RealSpace{2,Float64}
        @test rs_sv.x == SVector(1.0, 2.0)

        rs_tp = RealSpace((1.0, 2.0, 3.0))
        @test rs_tp isa RealSpace{3,Float64}
        @test rs_tp.x == SVector(1.0, 2.0, 3.0)
    end

    @testset "LatticeCoord construction" begin
        # default sublattice = 1
        lc1 = LatticeCoord((2, 3))
        @test lc1 isa LatticeCoord{2}
        @test lc1.cell == (2, 3)
        @test lc1.sublattice == 1

        # explicit sublattice
        lc2 = LatticeCoord((2, 3), 2)
        @test lc2.sublattice == 2

        # 1D
        lc_1d = LatticeCoord((5,))
        @test lc_1d isa LatticeCoord{1}
        @test lc_1d.cell == (5,)
    end

    @testset "HigherDimCoord construction" begin
        h = HigherDimCoord{2}(SVector(1.0, 2.0, 3.0, 4.0, 5.0))
        @test h isa HigherDimCoord{2,5,Float64}
        @test h.hyper == SVector(1.0, 2.0, 3.0, 4.0, 5.0)
    end

    @testset "identity conversions (lat-agnostic)" begin
        rs = RealSpace((1.0, 2.0))
        lc = LatticeCoord((3, 4))
        h = HigherDimCoord{2}(SVector(1.0, 2.0, 3.0, 4.0, 5.0))

        @test to_real(nothing, rs) === rs
        @test to_lattice(nothing, lc) === lc
        @test to_hyper(nothing, h) === h
    end
end
