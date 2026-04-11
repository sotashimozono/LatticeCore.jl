using LatticeCore
using Plots
using StaticArrays
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

@testset "diffraction_pattern on BraggPeakSet" begin
    # Build a tiny synthetic 1D BraggPeakSet (3 peaks: Γ + 2 side
    # peaks). We can do this by hand — the physics has already been
    # exercised in the QuasiCrystal test suite; here we just want
    # the plot method to return a `Plots.Plot`.
    peaks_1d = [SVector(0.0), SVector(1.0), SVector(-1.0)]
    intensities_1d = [1.0, 0.3, 0.3]
    hyper_indices_1d = [(0, 0), (1, 0), (-1, 0)]
    bps1 = BraggPeakSet{1,2,Float64}(peaks_1d, intensities_1d, hyper_indices_1d)

    p1 = diffraction_pattern(bps1; title="1D test")
    @test p1 isa Plots.Plot

    p1_log = diffraction_pattern(bps1; log_intensity=true)
    @test p1_log isa Plots.Plot

    # 2D case — 5 peaks in a plus + Γ configuration.
    peaks_2d = [
        SVector(0.0, 0.0),
        SVector(1.0, 0.0),
        SVector(-1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.0, -1.0),
    ]
    intensities_2d = [1.0, 0.5, 0.5, 0.5, 0.5]
    hyper_indices_2d = [
        (0, 0, 0, 0), (1, 0, 0, 0), (-1, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0)
    ]
    bps2 = BraggPeakSet{2,4,Float64}(peaks_2d, intensities_2d, hyper_indices_2d)

    p2 = diffraction_pattern(bps2; title="2D test")
    @test p2 isa Plots.Plot

    p2_log = diffraction_pattern(bps2; log_intensity=true, marker_scale=15.0)
    @test p2_log isa Plots.Plot

    # Uniform-intensity degenerate case: single peak at Γ.
    bps_single = BraggPeakSet{2,4,Float64}([SVector(0.0, 0.0)], [1.0], [(0, 0, 0, 0)])
    @test diffraction_pattern(bps_single) isa Plots.Plot
end
