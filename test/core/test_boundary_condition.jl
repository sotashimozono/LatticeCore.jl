using LatticeCore
using Test

@testset "BoundaryCondition type hierarchy" begin
    @test AbstractAxisBC <: Any
    @test PeriodicAxis <: AbstractAxisBC
    @test OpenAxis <: AbstractAxisBC
    @test TwistedAxis <: AbstractAxisBC

    @test AbstractBoundaryModifier <: Any
    @test NoModifier <: AbstractBoundaryModifier
    @test SSD <: AbstractBoundaryModifier

    @test LatticeBoundary <: AbstractBoundaryCondition
end

@testset "LatticeBoundary construction" begin
    b1 = LatticeBoundary((PeriodicAxis(),))
    @test b1 isa LatticeBoundary{1}
    @test length(b1.axes) == 1
    @test b1.axes[1] isa PeriodicAxis
    @test b1.modifier isa NoModifier

    b2 = LatticeBoundary((PeriodicAxis(), OpenAxis()))
    @test b2 isa LatticeBoundary{2}
    @test b2.axes[1] isa PeriodicAxis
    @test b2.axes[2] isa OpenAxis

    b3 = LatticeBoundary((PeriodicAxis(), PeriodicAxis()), SSD(8.0))
    @test b3.modifier isa SSD
    @test b3.modifier.L == 8.0
end

@testset "apply_axis_bc" begin
    @testset "PeriodicAxis wraps via mod1" begin
        @test apply_axis_bc(PeriodicAxis(), 0, 5) == (5, true)
        @test apply_axis_bc(PeriodicAxis(), 1, 5) == (1, true)
        @test apply_axis_bc(PeriodicAxis(), 5, 5) == (5, true)
        @test apply_axis_bc(PeriodicAxis(), 6, 5) == (1, true)
        @test apply_axis_bc(PeriodicAxis(), -1, 5) == (4, true)
    end

    @testset "OpenAxis rejects out-of-range" begin
        @test apply_axis_bc(OpenAxis(), 0, 5) == (0, false)
        @test apply_axis_bc(OpenAxis(), 1, 5) == (1, true)
        @test apply_axis_bc(OpenAxis(), 3, 5) == (3, true)
        @test apply_axis_bc(OpenAxis(), 5, 5) == (5, true)
        @test apply_axis_bc(OpenAxis(), 6, 5) == (6, false)
    end

    @testset "TwistedAxis wraps like PBC" begin
        ax = TwistedAxis(π / 4)
        @test apply_axis_bc(ax, 6, 5) == (1, true)
        @test apply_axis_bc(ax, 0, 5) == (5, true)
    end
end

@testset "axis_phase" begin
    # Non-twisted axes return unit phase.
    @test axis_phase(PeriodicAxis(), 6, 5) == complex(1.0, 0.0)
    @test axis_phase(OpenAxis(), 6, 5) == complex(1.0, 0.0)

    # TwistedAxis attaches exp(±iθ) when the raw index exits the range.
    θ = 0.3
    ax = TwistedAxis(θ)
    @test axis_phase(ax, 6, 5) ≈ cis(θ)
    @test axis_phase(ax, 0, 5) ≈ cis(-θ)
    @test axis_phase(ax, 3, 5) == complex(1.0, 0.0)   # interior: no phase
end

@testset "bond_weight" begin
    # NoModifier always returns 1.0; no lattice actually needed.
    @test bond_weight(NoModifier(), nothing, 1, 2) == 1.0
end
