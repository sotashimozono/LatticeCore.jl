module LatticeCorePlotsExt

using LatticeCore
using LinearAlgebra
using Plots
using StaticArrays

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

# ---- diffraction_pattern --------------------------------------------

function _peak_arrays(ml::LatticeCore.AbstractMomentumLattice{DPhys,T}) where {DPhys,T}
    N = num_k_points(ml)
    ks = [k_point(ml, i) for i in 1:N]
    Is = if ml isa LatticeCore.BraggPeakSet
        copy(ml.intensities)
    else
        ones(T, N)
    end
    return ks, Is
end

function _intensity_for_plot(Is::AbstractVector{T}, log_intensity::Bool) where {T}
    if !log_intensity
        return Float64.(Is)
    end
    Imax = maximum(Is)
    Imax <= 0 && return fill(-Inf, length(Is))
    floor = 1e-6
    return [I > 0 ? log10(max(I / Imax, floor)) : log10(floor) for I in Is]
end

# 1D diffraction pattern: stem plot of intensity vs k.
function LatticeCore.diffraction_pattern(
    ml::LatticeCore.AbstractMomentumLattice{1,T};
    title="Diffraction pattern",
    log_intensity::Bool=false,
    color=:steelblue,
    lw=1.5,
    ms=4,
    kwargs...,
) where {T}
    ks, Is = _peak_arrays(ml)
    xs = [Float64(k[1]) for k in ks]
    ys = _intensity_for_plot(Is, log_intensity)

    p = Plots.plot(;
        xlabel="k",
        ylabel=log_intensity ? "log₁₀(I / I_max)" : "I(k)",
        title=title,
        legend=false,
        kwargs...,
    )
    seg_x = Float64[]
    seg_y = Float64[]
    base = log_intensity ? minimum(ys) : 0.0
    for (x, y) in zip(xs, ys)
        push!(seg_x, x, x, NaN)
        push!(seg_y, base, y, NaN)
    end
    Plots.plot!(p, seg_x, seg_y; color=color, lw=lw, label="")
    Plots.scatter!(p, xs, ys; mc=color, ms=ms, markerstrokewidth=0, label="")
    return p
end

# 2D diffraction pattern: scatter of (kx, ky) with marker size
# proportional to intensity.
function LatticeCore.diffraction_pattern(
    ml::LatticeCore.AbstractMomentumLattice{2,T};
    title="Diffraction pattern",
    log_intensity::Bool=false,
    marker_scale::Real=12.0,
    color=:viridis,
    kwargs...,
) where {T}
    ks, Is = _peak_arrays(ml)
    xs = [Float64(k[1]) for k in ks]
    ys = [Float64(k[2]) for k in ks]
    zs = _intensity_for_plot(Is, log_intensity)

    # Map intensities to marker sizes (raw or log).
    zmax = maximum(zs)
    zmin = minimum(zs)
    rng = zmax - zmin
    ms = if rng > 0
        [1.0 + marker_scale * (z - zmin) / rng for z in zs]
    else
        fill(Float64(marker_scale) / 2, length(zs))
    end

    p = Plots.scatter(
        xs,
        ys;
        ms=ms,
        marker_z=zs,
        color=color,
        markerstrokewidth=0,
        aspect_ratio=:equal,
        xlabel="kₓ",
        ylabel="k_y",
        title=title,
        label="",
        colorbar_title=log_intensity ? "log₁₀(I / I_max)" : "I(k)",
        kwargs...,
    )
    return p
end

# ---- diffraction_pattern from (lat, state, q_grid) ------------------
#
# Issue #27: the form `diffraction_pattern(lat, state, q_grid)` runs
# `structure_factor` on the supplied grid (so the FFT / NUFFT fast
# paths kick in transparently when the matching extension is loaded)
# and renders the result. For 2D regular meshes we draw a heatmap on
# the natural rectangular grid; otherwise we delegate to the existing
# scatter / stem methods.

function _select_intensity(Sks, intensity::Symbol)
    if intensity === :S
        return Sks
    elseif intensity === Symbol("|S|²") || intensity === :Smag2 || intensity === :intensity
        return Sks  # `structure_factor` already returns |S|² / N
    else
        throw(ArgumentError("intensity must be one of (:|S|², :S); got $(intensity)"))
    end
end

# 1D → stem plot of S(k) vs k.
function LatticeCore.diffraction_pattern(
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice{1,T};
    title="Diffraction pattern",
    intensity::Symbol=Symbol("|S|²"),
    log_intensity::Bool=false,
    color=:steelblue,
    lw=1.5,
    ms=4,
    kwargs...,
) where {T}
    Sks = _select_intensity(structure_factor(lat, state, ml), intensity)
    M = num_k_points(ml)
    xs = [Float64(k_point(ml, i)[1]) for i in 1:M]
    ys = _intensity_for_plot(Sks, log_intensity)

    p = Plots.plot(;
        xlabel="k",
        ylabel=log_intensity ? "log₁₀(I / I_max)" : "I(k)",
        title=title,
        legend=false,
        kwargs...,
    )
    seg_x = Float64[]
    seg_y = Float64[]
    base = log_intensity ? minimum(ys) : 0.0
    perm = sortperm(xs)
    for j in perm
        push!(seg_x, xs[j], xs[j], NaN)
        push!(seg_y, base, ys[j], NaN)
    end
    Plots.plot!(p, seg_x, seg_y; color=color, lw=lw, label="")
    Plots.scatter!(p, xs[perm], ys[perm]; mc=color, ms=ms, markerstrokewidth=0, label="")
    return p
end

# 2D regular mesh → heatmap. We use the mesh dims to reshape the flat
# `structure_factor` output back into a (Nx, Ny) grid; the kx / ky
# axes are the unique values on each axis.
function LatticeCore.diffraction_pattern(
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector,
    ml::LatticeCore.PeriodicMomentumLattice{2,T};
    title="Diffraction pattern",
    intensity::Symbol=Symbol("|S|²"),
    log_intensity::Bool=false,
    color=:viridis,
    kwargs...,
) where {T}
    Sks = _select_intensity(structure_factor(lat, state, ml), intensity)
    Nx, Ny = ml.mesh

    # `gamma_centered` and `monkhorst_pack` enumerate k-points via
    # `CartesianIndices(mesh)` in column-major order, so reshape lines
    # the values up with (kx, ky) on a Cartesian grid.
    grid = reshape(_intensity_for_plot(Sks, log_intensity), Nx, Ny)
    kxs = [Float64(k_point(ml, i)[1]) for i in 1:Nx]
    kys = [Float64(k_point(ml, (j - 1) * Nx + 1)[2]) for j in 1:Ny]

    return Plots.heatmap(
        kxs,
        kys,
        permutedims(grid);  # heatmap wants (Ny, Nx) → rows are y
        aspect_ratio=:equal,
        xlabel="kₓ",
        ylabel="k_y",
        title=title,
        color=color,
        colorbar_title=log_intensity ? "log₁₀(I / I_max)" : "I(k)",
        kwargs...,
    )
end

# 2D irregular momentum lattice (e.g. BraggPeakSet) → reuse the scatter
# method on a freshly-built BraggPeakSet of computed intensities.
function LatticeCore.diffraction_pattern(
    lat::LatticeCore.AbstractLattice,
    state::AbstractVector,
    ml::LatticeCore.AbstractMomentumLattice{2,T};
    title="Diffraction pattern",
    intensity::Symbol=Symbol("|S|²"),
    log_intensity::Bool=false,
    marker_scale::Real=12.0,
    color=:viridis,
    kwargs...,
) where {T}
    Sks = Float64.(_select_intensity(structure_factor(lat, state, ml), intensity))
    M = num_k_points(ml)
    ks = [k_point(ml, i) for i in 1:M]
    bps = LatticeCore.BraggPeakSet{2,2,Float64}(
        [SVector{2,Float64}(Float64(k[1]), Float64(k[2])) for k in ks],
        Sks,
        [(0, 0) for _ in 1:M],
    )
    return LatticeCore.diffraction_pattern(
        bps;
        title=title,
        log_intensity=log_intensity,
        marker_scale=marker_scale,
        color=color,
        kwargs...,
    )
end

end # module LatticeCorePlotsExt
