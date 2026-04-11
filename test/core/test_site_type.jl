using LatticeCore
using LinearAlgebra
using StaticArrays
using Random
using Test

@testset "AbstractSiteType" begin
    @testset "type hierarchy" begin
        @test IsingSite <: AbstractSiteType
        @test PottsSite <: AbstractSiteType
        @test XYSite <: AbstractSiteType
        @test HeisenbergSite <: AbstractSiteType
        @test EmptySite <: AbstractSiteType
    end

    @testset "element_type default is VertexCenter" begin
        @test element_type(IsingSite()) isa VertexCenter
        @test element_type(XYSite()) isa VertexCenter
        @test element_type(HeisenbergSite()) isa VertexCenter
        @test element_type(PottsSite(3)) isa VertexCenter
        @test element_type(EmptySite()) isa VertexCenter
    end

    @testset "IsingSite" begin
        s = IsingSite()
        @test s isa IsingSite{Int8}
        @test state_type(s) === Int8
        @test zero_state(s) === Int8(0)
        @test domain(s) == (Int8(-1), Int8(1))

        rng = MersenneTwister(42)
        for _ in 1:20
            v = random_state(rng, s)
            @test v isa Int8
            @test v == Int8(-1) || v == Int8(1)
        end

        # Custom storage type
        s_int = IsingSite{Int}()
        @test state_type(s_int) === Int
        @test random_state(MersenneTwister(0), s_int) isa Int
    end

    @testset "PottsSite" begin
        s = PottsSite(3)
        @test s isa PottsSite{3,Int8}
        @test state_type(s) === Int8
        @test domain(s) == (Int8(1), Int8(2), Int8(3))

        rng = MersenneTwister(42)
        for _ in 1:30
            v = random_state(rng, s)
            @test v isa Int8
            @test 1 <= v <= 3
        end
    end

    @testset "XYSite" begin
        s = XYSite()
        @test s isa XYSite{Float64}
        @test state_type(s) === Float64
        @test zero_state(s) === 0.0

        rng = MersenneTwister(42)
        for _ in 1:20
            θ = random_state(rng, s)
            @test θ isa Float64
            @test 0.0 <= θ < 2π
        end
    end

    @testset "HeisenbergSite (unit vectors)" begin
        s = HeisenbergSite()
        @test s isa HeisenbergSite{Float64}
        @test state_type(s) === SVector{3,Float64}

        rng = MersenneTwister(42)
        for _ in 1:30
            v = random_state(rng, s)
            @test v isa SVector{3,Float64}
            @test isapprox(norm(v), 1.0; atol=1e-12)
        end
    end

    @testset "EmptySite" begin
        s = EmptySite()
        @test state_type(s) === Nothing
        @test random_state(MersenneTwister(0), s) === nothing
        @test zero_state(s) === nothing
    end
end
