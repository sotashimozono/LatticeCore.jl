using LatticeCore
using Test

# ---- Minimal mock of an infinite lattice, local to this test file ----
#
# A proper infinite lattice like `InfiniteFibonacci` is out of scope
# for LatticeCore itself; it will live in a downstream package
# (QuasiCrystal.jl / a future LSystemAdapter). We only need a
# *contract check* here: `size_trait = InfiniteSize()` must make
# `is_finite` false and trip `require_finite`.

struct _InfiniteChain <: AbstractLattice{1,Float64} end

LatticeCore.size_trait(::_InfiniteChain) = InfiniteSize()
LatticeCore.boundary(::_InfiniteChain) = LatticeBoundary((PeriodicAxis(),), NoModifier())
LatticeCore.site_layout(::_InfiniteChain) = UniformLayout(IsingSite())

# And a minimal "infinite abstract" + matching `materialize` overload,
# demonstrating the intended usage pattern from the 06 chapter.

struct _InfiniteChainAbstract end

LatticeCore.materialize(::_InfiniteChainAbstract; depth::Int) =
    LineLattice(depth, PeriodicAxis())

@testset "LazyInfinite" begin
    @testset "require_finite on finite lattices is a no-op" begin
        @test require_finite(LineLattice(5)) === nothing
        @test require_finite(SimpleSquareLattice(3, 3)) === nothing
        @test require_finite(LineLattice(5, OpenAxis())) === nothing
    end

    @testset "require_finite throws on InfiniteSize" begin
        inf_lat = _InfiniteChain()
        @test is_finite(inf_lat) == false
        @test_throws ArgumentError require_finite(inf_lat)
    end

    @testset "materialize turns an infinite abstract into a finite lattice" begin
        abs_chain = _InfiniteChainAbstract()

        lat = materialize(abs_chain; depth=10)
        @test lat isa LineLattice
        @test num_sites(lat) == 10
        @test is_finite(lat) == true

        # The materialised lattice passes the MC guard.
        @test require_finite(lat) === nothing
    end

    @testset "materialize is a generic function with no default" begin
        # Calling it on a bare type with no specialised method should
        # raise a MethodError, because `materialize` deliberately has
        # no fallback implementation.
        @test_throws MethodError materialize(Ref(nothing); depth=1)
    end
end
