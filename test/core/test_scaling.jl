using LatticeCore
using Test

# A minimal lattice that opts out of scaling (the default), and one that opts in,
# so the generic fallbacks and `size_sequence` can both be exercised without
# depending on a concrete lattice package.
struct _NoScaleLattice <: AbstractLattice{2,Float64} end

struct _StepLattice <: AbstractLattice{2,Float64}
    n::Int
end
LatticeCore.scaling_rule(::_StepLattice) = LinearScaling(2)
LatticeCore.rescale(l::_StepLattice, k::Integer=1) = _StepLattice(l.n * 2^k)
function LatticeCore.cell_partition(l::_StepLattice, k::Integer=1)
    return [collect(((c - 1) * 2 ^ k + 1):(c * 2 ^ k)) for c in 1:(l.n ÷ 2 ^ k)]
end

@testset "scaling rule traits" begin
    @test NoScaling() isa AbstractScalingRule
    @test LinearScaling(2) isa AbstractScalingRule
    @test SubstitutionScaling() isa AbstractScalingRule

    @test LinearScaling(3).factor == 3
    @test SubstitutionScaling().depth_step == 1
    @test SubstitutionScaling(2).depth_step == 2

    # a scale step must actually change the scale
    @test_throws ArgumentError LinearScaling(1)
    @test_throws ArgumentError LinearScaling(0)
    @test_throws ArgumentError SubstitutionScaling(0)
    @test_throws ArgumentError SubstitutionScaling(-1)
end

@testset "defaults: a lattice opts out unless it says otherwise" begin
    l = _NoScaleLattice()
    @test scaling_rule(l) === NoScaling()
    # no silent no-op: an unsupported scale change is an error, not the identity
    @test_throws ArgumentError rescale(l)
    @test_throws ArgumentError rescale(l, 3)
    @test_throws ArgumentError cell_partition(l)
end

@testset "size_sequence rides rescale" begin
    l = _StepLattice(1)
    @test scaling_rule(l) === LinearScaling(2)

    @test rescale(l, 0).n == 1          # k = 0 is the lattice itself
    @test rescale(l).n == 2             # default step is one
    @test rescale(l, 3).n == 8

    seq = size_sequence(l, 4)
    @test length(seq) == 5
    @test [s.n for s in seq] == [1, 2, 4, 8, 16]
end

@testset "cell_partition partitions the sites" begin
    l = _StepLattice(8)
    for k in 1:3
        groups = cell_partition(l, k)
        flat = reduce(vcat, groups)
        @test sort(flat) == collect(1:8)          # a partition: covers, no repeats
        @test length(groups) == 8 ÷ 2^k
        @test all(g -> length(g) == 2^k, groups)
    end
end
