module LatticeCorePlotsExt

using LatticeCore
using LinearAlgebra
using Plots

# ---- Common helpers --------------------------------------------------

"""
    _site_xy(lat) → (xs, ys)

Collect the x and y coordinates of every site as a pair of
`Vector{Float64}`. 1D lattices get `ys = zeros(N)`.
"""
function _site_xy(lat::LatticeCore.AbstractLattice{1})
    N = num_sites(lat)
    xs = [Float64(position(lat, i)[1]) for i in 1:N]
    ys = zeros(N)
    return xs, ys
end

function _site_xy(lat::LatticeCore.AbstractLattice{2})
    N = num_sites(lat)
    xs = [Float64(position(lat, i)[1]) for i in 1:N]
    ys = [Float64(position(lat, i)[2]) for i in 1:N]
    return xs, ys
end

"""
    _bond_segments(lat) → (seg_x, seg_y)

Collect bond line segments as flat `Vector{Float64}` with `NaN`
separators — the format `Plots.plot!` expects for drawing a
disjoint union of line segments in a single call.

Each bond's endpoints are computed as `src + b.vector` rather than
`position(lat, b.j)`, which lets the plot respect the wrapped
displacement vector stored on each `Bond` and avoids "long lines"
artefacts on periodic lattices.
"""
function _bond_segments(lat::LatticeCore.AbstractLattice{1})
    seg_x = Float64[]
    seg_y = Float64[]
    for b in bonds(lat)
        src = position(lat, b.i)
        dst_local = src + b.vector
        push!(seg_x, Float64(src[1]), Float64(dst_local[1]), NaN)
        push!(seg_y, 0.0, 0.0, NaN)
    end
    return seg_x, seg_y
end

function _bond_segments(lat::LatticeCore.AbstractLattice{2})
    seg_x = Float64[]
    seg_y = Float64[]
    for b in bonds(lat)
        src = position(lat, b.i)
        dst_local = src + b.vector
        push!(seg_x, Float64(src[1]), Float64(dst_local[1]), NaN)
        push!(seg_y, Float64(src[2]), Float64(dst_local[2]), NaN)
    end
    return seg_x, seg_y
end

# ---- plot_sites! -----------------------------------------------------

function LatticeCore.plot_sites!(
    p, lat::LatticeCore.AbstractLattice{D}; marker_size=4, color=:auto
) where {D}
    xs, ys = _site_xy(lat)
    N = length(xs)

    if color === :auto && num_sublattices(lat) > 1
        sub_ids = [sublattice(lat, i) for i in 1:N]
        Plots.scatter!(p, xs, ys; ms=marker_size, group=sub_ids, markerstrokewidth=0)
    else
        mc = color === :auto ? :steelblue : color
        Plots.scatter!(p, xs, ys; ms=marker_size, mc=mc, markerstrokewidth=0, label="")
    end
    return p
end

# ---- plot_bonds! -----------------------------------------------------

function LatticeCore.plot_bonds!(
    p, lat::LatticeCore.AbstractLattice{D}; color=:black, lw=1.0
) where {D}
    seg_x, seg_y = _bond_segments(lat)
    return Plots.plot!(p, seg_x, seg_y; color=color, lw=lw, label="")
end

# ---- plot_lattice (1D) ----------------------------------------------

function LatticeCore.plot_lattice(
    lat::LatticeCore.AbstractLattice{1,T};
    show_bonds::Bool=true,
    show_sites::Bool=true,
    site_size=6,
    site_color=:auto,
    bond_color=:black,
    bond_width=1.5,
    title=nothing,
    kwargs...,
) where {T}
    p = Plots.plot(;
        xlabel="x", ylims=(-0.5, 0.5), yticks=[], legend=:topright, title=title, kwargs...
    )
    show_bonds && plot_bonds!(p, lat; color=bond_color, lw=bond_width)
    show_sites && plot_sites!(p, lat; marker_size=site_size, color=site_color)
    return p
end

# ---- plot_lattice (2D) ----------------------------------------------

function LatticeCore.plot_lattice(
    lat::LatticeCore.AbstractLattice{2,T};
    show_bonds::Bool=true,
    show_sites::Bool=true,
    site_size=4,
    site_color=:auto,
    bond_color=:black,
    bond_width=1.0,
    aspect_ratio=:equal,
    title=nothing,
    kwargs...,
) where {T}
    p = Plots.plot(;
        aspect_ratio=aspect_ratio,
        xlabel="x",
        ylabel="y",
        legend=:topright,
        title=title,
        kwargs...,
    )
    show_bonds && plot_bonds!(p, lat; color=bond_color, lw=bond_width)
    show_sites && plot_sites!(p, lat; marker_size=site_size, color=site_color)
    return p
end

end # module LatticeCorePlotsExt
