using LatticeCore
using Plots
using Test

# Loading Plots triggers the LatticeCorePlotsExt extension, which
# installs methods on `LatticeCore.plot_lattice`, `plot_bonds!`, and
# `plot_sites!`. The tests below verify that those methods are
# available and that they return a `Plots.Plot` on every shipped
# reference lattice.

const _LAT_CORE_EXT = Base.get_extension(LatticeCore, :LatticeCorePlotsExt)

@testset "LatticeCorePlotsExt loads when Plots is imported" begin
    @test _LAT_CORE_EXT !== nothing

    # `plot_lattice` should now resolve to at least one concrete method
    # for the 1D and 2D reference lattices.
    @test !isempty(methods(plot_lattice, (AbstractLattice{1,Float64},)))
    @test !isempty(methods(plot_lattice, (AbstractLattice{2,Float64},)))
end

@testset "plot_lattice on LineLattice" begin
    lat = LineLattice(6, PeriodicAxis())
    p = plot_lattice(lat; title="LineLattice PBC 6")
    @test p isa Plots.Plot

    # Without bonds
    p_no_bonds = plot_lattice(lat; show_bonds=false)
    @test p_no_bonds isa Plots.Plot

    # Custom site colour
    p_colored = plot_lattice(lat; site_color=:red)
    @test p_colored isa Plots.Plot
end

@testset "plot_lattice on SimpleSquareLattice" begin
    lat = SimpleSquareLattice(3, 3, PeriodicAxis())
    p = plot_lattice(lat; title="Square PBC 3x3")
    @test p isa Plots.Plot

    # Cylinder: mixed BC still works
    cyl = SimpleSquareLattice(3, 4, LatticeBoundary((PeriodicAxis(), OpenAxis())))
    p_cyl = plot_lattice(cyl; title="Cylinder 3x4")
    @test p_cyl isa Plots.Plot

    # show_sites = false still produces a plot
    p_bonds_only = plot_lattice(lat; show_sites=false)
    @test p_bonds_only isa Plots.Plot
end

@testset "plot_bonds! / plot_sites! on an existing plot" begin
    lat = SimpleSquareLattice(4, 4, PeriodicAxis())

    p = Plots.plot(; aspect_ratio=:equal, label="")
    plot_sites!(p, lat; marker_size=5, color=:steelblue)
    plot_bonds!(p, lat; color=:gray, lw=0.8)
    @test p isa Plots.Plot
end

@testset "wrapped bond vectors keep bond lines short under PBC" begin
    # On a 3x3 PBC square the wrapped bond at the corner is (−1, 0),
    # not (2, 0). _bond_segments must respect that, so the total length
    # of drawn segments is bounded by (unit) × (number of bonds × 2
    # endpoints) rather than exploding across the whole sample.
    lat = SimpleSquareLattice(3, 3, PeriodicAxis())
    ext = _LAT_CORE_EXT
    seg_x, seg_y = ext._bond_segments(lat)

    # Each segment is (src, dst, NaN) = 3 entries.
    @test length(seg_x) == 3 * length(collect(bonds(lat)))
    @test length(seg_y) == 3 * length(collect(bonds(lat)))

    # Every drawn segment has unit length in one of (x, y).
    for i in 1:3:length(seg_x)
        dx = seg_x[i + 1] - seg_x[i]
        dy = seg_y[i + 1] - seg_y[i]
        @test isapprox(hypot(dx, dy), 1.0; atol=1e-10)
    end
end
