using CairoMakie          # loads Makie -> triggers LatticeCoreMakieExt
using Test

const MAKIE_EXT = Base.get_extension(LatticeCore, :LatticeCoreMakieExt)

renders(fig) = mktemp() do path, io
    close(io)
    p = path * ".png"
    CairoMakie.save(p, fig)
    ok = isfile(p) && filesize(p) > 0
    rm(p; force=true)
    return ok
end

@testset "LatticeCoreMakieExt" begin
    @testset "extension loads with Makie" begin
        @test MAKIE_EXT !== nothing
    end

    @testset "plot_lattice backend dispatch: 2D and 1D, renderable" begin
        for lat in (SimpleSquareLattice(4, 4), LineLattice(8, OpenAxis()))
            fig = plot_lattice(lat; backend=MakieBackend())
            @test fig isa Makie.Figure
            @test renders(fig)
        end
        sq = SimpleSquareLattice(4, 4)
        @test plot_lattice(sq; backend=:makie, highlight_bonds=[1, 2, 3]) isa Makie.Figure
        @test plot_lattice(sq; backend=MakieBackend(), colorby=:none) isa Makie.Figure
    end

    @testset "backend registry & symbol mapping" begin
        @test MakieBackend() in LatticeCore._LOADED_PLOT_BACKENDS   # registered on load
        @test default_plot_backend() isa AbstractPlotBackend
        @test LatticeCore._as_backend(:makie) === MakieBackend()
        @test LatticeCore._as_backend(:plots) === PlotsBackend()
        @test LatticeCore._as_backend(MakieBackend()) === MakieBackend()
        @test_throws ArgumentError plot_lattice(SimpleSquareLattice(3, 3); backend=:nope)
    end

    @testset "makie_state: renderable, arrows, validation" begin
        lat = SimpleSquareLattice(5, 5)
        N = num_sites(lat)
        fig = makie_state(lat, Float64.(1:N))
        @test fig isa Makie.Figure
        @test renders(fig)
        @test makie_state(lat, [2π * i / N for i in 1:N]; arrows=true) isa Makie.Figure
        @test_throws DimensionMismatch makie_state(lat, Float64.(1:(N - 1)))
    end

    @testset "makie_structure_factor: figure + S(k) numerics" begin
        lat = SimpleSquareLattice(6, 6)
        N = num_sites(lat)
        @test makie_structure_factor(lat, ones(N); resolution=16) isa Makie.Figure

        R = 21
        ks, S = MAKIE_EXT._structure_factor_grid(
            lat, ones(N); k_range=(-π, π), resolution=R
        )
        c = (R + 1) ÷ 2
        @test abs(ks[c]) < 1e-12                       # centre is k = 0
        @test S[c, c] ≈ float(N) atol = 1e-8           # uniform state ⇒ Γ peak = N
        @test S ≈ reverse(S) rtol = 1e-10              # S(k) = S(-k)
        @test all(S .≥ -1e-12)
        @test_throws DimensionMismatch makie_structure_factor(lat, ones(N - 1))
        @test_throws ArgumentError makie_structure_factor(lat, ones(N); resolution=0)
    end
end
