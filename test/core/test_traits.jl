using LatticeCore
using Test

@testset "Traits" begin
    @testset "TopologyTrait" begin
        @test TopologyTrait{:Square}() isa TopologyTrait{:Square}
        @test TopologyTrait{:Honeycomb} <: Any
    end

    @testset "Periodicity" begin
        @test Periodic() isa Periodic
        @test Aperiodic() isa Aperiodic
    end

    @testset "AbstractReciprocalSupport" begin
        @test HasReciprocal <: AbstractReciprocalSupport
        @test HasFourierModule <: AbstractReciprocalSupport
        @test NoReciprocal <: AbstractReciprocalSupport
        @test HasReciprocal() isa AbstractReciprocalSupport
        @test HasFourierModule() isa AbstractReciprocalSupport
        @test NoReciprocal() isa AbstractReciprocalSupport
    end

    @testset "AbstractSizeTrait" begin
        @test FiniteSize <: AbstractSizeTrait
        @test InfiniteSize <: AbstractSizeTrait
        @test QuasiInfiniteSize <: AbstractSizeTrait

        finite2d = FiniteSize((4, 5))
        @test finite2d isa FiniteSize{2}
        @test finite2d.dims == (4, 5)

        @test InfiniteSize() isa InfiniteSize

        qi = QuasiInfiniteSize(3.5)
        @test qi isa QuasiInfiniteSize{Float64}
        @test qi.cutoff == 3.5
    end
end
