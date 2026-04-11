using LatticeCore
using Test

@testset "AbstractSiteLayout" begin
    @testset "type hierarchy" begin
        @test UniformLayout <: AbstractSiteLayout
        @test SublatticeLayout <: AbstractSiteLayout
        @test ExplicitLayout <: AbstractSiteLayout
    end

    @testset "UniformLayout" begin
        l = UniformLayout(IsingSite())
        @test l isa UniformLayout{IsingSite{Int8}}
        @test l.site_type isa IsingSite

        # Every site returns the same site type regardless of i
        @test site_type(l, 1) === l.site_type
        @test site_type(l, 42) === l.site_type
        @test site_type(l, 10_000) === l.site_type
    end

    @testset "SublatticeLayout (mixed-spin)" begin
        # Honeycomb-style mixed layout: A = Ising, B = XY
        sublattice_of = [1, 2, 1, 2, 1, 2]  # ABABAB
        l = SublatticeLayout((IsingSite(), XYSite()), sublattice_of)
        @test l isa SublatticeLayout{2}

        @test site_type(l, 1) isa IsingSite   # A
        @test site_type(l, 2) isa XYSite      # B
        @test site_type(l, 5) isa IsingSite
        @test site_type(l, 6) isa XYSite
    end

    @testset "ExplicitLayout" begin
        types = AbstractSiteType[IsingSite(), IsingSite(), XYSite(), HeisenbergSite()]
        l = ExplicitLayout(types)
        @test site_type(l, 1) isa IsingSite
        @test site_type(l, 3) isa XYSite
        @test site_type(l, 4) isa HeisenbergSite
    end
end
