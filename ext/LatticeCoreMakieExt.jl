module LatticeCoreMakieExt

using LatticeCore
using Makie

"""
    LatticeCoreMakieExt

Makie backend for LatticeCore, loaded automatically once `Makie` is in scope
(e.g. `using CairoMakie` / `GLMakie`). It is the Makie counterpart of
`LatticeCorePlotsExt`, implementing the generic `makie_*` entry points on any
`AbstractLattice`:

- `plot_lattice(lat; backend = MakieBackend())` — sites + bonds, sublattice
  colouring, bond highlight (the Makie method of the shared, backend-dispatched
  `plot_lattice`);
- [`makie_state`](@ref) — per-site scalar field as a colour-mapped scatter,
  optional in-plane spin arrows;
- [`makie_structure_factor`](@ref) — static `S(k)` heatmap (2D lattices).

Each returns a `Makie.Figure`.
"""
LatticeCoreMakieExt

# Register this backend so `plot_lattice(lat; backend = default_plot_backend())`
# resolves to Makie once the extension is loaded.
__init__() = LatticeCore._register_plot_backend(LatticeCore.MakieBackend())

# Real-space (x, y) of site i, padding 1D lattices to the y = 0 line.
@inline function _xy(lat, i)
    p = LatticeCore.position(lat, i)
    x = Float64(p[1])
    y = length(p) ≥ 2 ? Float64(p[2]) : 0.0
    return x, y
end

function _site_coords(lat)
    N = LatticeCore.num_sites(lat)
    xs = Vector{Float64}(undef, N)
    ys = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        xs[i], ys[i] = _xy(lat, i)
    end
    return xs, ys
end

# Bond segments drawn from site i along the per-bond wrapped displacement, so
# periodic samples show no boundary-crossing long lines.
function _bond_segments(lat, bs, indices)
    seg = Point2f[]
    for k in indices
        b = bs[k]
        xi, yi = _xy(lat, b.i)
        dx = Float64(b.vector[1])
        dy = length(b.vector) ≥ 2 ? Float64(b.vector[2]) : 0.0
        push!(seg, Point2f(xi, yi))
        push!(seg, Point2f(xi + dx, yi + dy))
    end
    return seg
end

function LatticeCore.plot_lattice(
    ::LatticeCore.MakieBackend,
    lat::LatticeCore.AbstractLattice;
    colorby::Symbol=:sublattice,
    highlight_bonds=nothing,
    markersize::Real=12,
    show_sites::Bool=true,
    figure=(;),
    axis=(;),
)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), axis...)

    bs = collect(LatticeCore.bonds(lat))
    linesegments!(
        ax, _bond_segments(lat, bs, eachindex(bs)); color=(:gray, 0.6), linewidth=1.0
    )
    if highlight_bonds !== nothing
        linesegments!(
            ax, _bond_segments(lat, bs, highlight_bonds); color=:red, linewidth=2.5
        )
    end

    if show_sites
        xs, ys = _site_coords(lat)
        if colorby === :sublattice && LatticeCore.num_sublattices(lat) > 1
            cols = [LatticeCore.sublattice(lat, i) for i in eachindex(xs)]
            scatter!(ax, xs, ys; color=cols, colormap=:tab10, markersize=markersize)
        else
            scatter!(ax, xs, ys; color=:steelblue, markersize=markersize)
        end
    end
    return fig
end

function LatticeCore.makie_state(
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector;
    colormap=:RdBu,
    arrows::Bool=false,
    markersize::Real=15,
    figure=(;),
    axis=(;),
)
    N = LatticeCore.num_sites(lat)
    length(state) == N || throw(
        DimensionMismatch("state has length $(length(state)) but lattice has $N sites")
    )
    xs, ys = _site_coords(lat)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), axis...)
    sc = scatter!(
        ax, xs, ys; color=Float64.(state), colormap=colormap, markersize=markersize
    )
    Colorbar(fig[1, 2], sc)
    if arrows
        seg = Point2f[]
        @inbounds for i in 1:N
            θ = Float64(state[i])
            dx, dy = 0.4 * cos(θ), 0.4 * sin(θ)
            push!(seg, Point2f(xs[i] - dx, ys[i] - dy))
            push!(seg, Point2f(xs[i] + dx, ys[i] + dy))
        end
        linesegments!(ax, seg; color=:black, linewidth=1.5)
    end
    return fig
end

# S(k) = |Σ_j state_j exp(-i k·r_j)|² / N on a resolution×resolution k-grid.
# Split out so the numerics can be tested without rendering. 2D lattices only.
function _structure_factor_grid(lat, state::AbstractVector; k_range, resolution::Int)
    N = LatticeCore.num_sites(lat)
    length(state) == N || throw(
        DimensionMismatch("state has length $(length(state)) but lattice has $N sites")
    )
    resolution ≥ 1 || throw(ArgumentError("resolution must be ≥ 1, got $resolution"))
    xs, ys = _site_coords(lat)
    st = ComplexF64.(state)
    ks = range(Float64(k_range[1]), Float64(k_range[2]); length=resolution)
    S = Matrix{Float64}(undef, resolution, resolution)
    @inbounds for a in 1:resolution
        kx = ks[a]
        for b in 1:resolution
            ky = ks[b]
            acc = zero(ComplexF64)
            for j in 1:N
                acc += st[j] * cis(-(kx * xs[j] + ky * ys[j]))
            end
            S[a, b] = abs2(acc) / N
        end
    end
    return ks, S
end

function LatticeCore.makie_structure_factor(
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector;
    k_range=(-π, π),
    resolution::Int=200,
    colormap=:viridis,
    figure=(;),
    axis=(;),
)
    ks, S = _structure_factor_grid(lat, state; k_range=k_range, resolution=resolution)
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="kₓ", ylabel="k_y", axis...)
    hm = heatmap!(ax, ks, ks, S; colormap=colormap)
    Colorbar(fig[1, 2], hm)
    return fig
end

end # module LatticeCoreMakieExt
