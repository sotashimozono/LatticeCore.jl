"""
    plot_lattice(lat::AbstractLattice{D, T}; kwargs...)

Generic visualisation entry point for any `AbstractLattice`. Returns
a `Plots.Plot`.

This function is a stub in the core module. Concrete methods live in
the [`LatticeCorePlotsExt`](@ref) package extension and are loaded
automatically when a client module imports `Plots`:

```julia
using Plots, LatticeCore
plot_lattice(SimpleSquareLattice(4, 4))
```

# Common keyword arguments

- `show_bonds::Bool = true` — draw bonds as line segments.
- `show_sites::Bool = true` — draw sites as scatter points.
- `site_size::Real` — marker size for scatter points.
- `site_color` — `:auto` (colour by sublattice when `num_sublattices >
  1`), or any `Plots.jl` colour value.
- `bond_color`, `bond_width` — bond styling.
- Any other keyword is forwarded to the underlying `Plots.plot` call
  (e.g. `title`, `xlabel`, `xlims`).

Calling `plot_lattice` without loading a backend (`Plots` or `CairoMakie`)
first raises an error that directs the user to the extension.
"""
function plot_lattice(lat::AbstractLattice; backend=default_plot_backend(), kwargs...)
    return plot_lattice(_as_backend(backend), lat; kwargs...)
end

# ---- Plotting backends ----------------------------------------------

"""
    AbstractPlotBackend

Base type for the plotting backends that the `plot_*` functions dispatch on.
Concrete singletons [`PlotsBackend`](@ref) and [`MakieBackend`](@ref) are
provided by `LatticeCorePlotsExt` / `LatticeCoreMakieExt` once `Plots` /
`Makie` is loaded.
"""
abstract type AbstractPlotBackend end

"""    PlotsBackend() — dispatch tag selecting the `Plots.jl` backend."""
struct PlotsBackend <: AbstractPlotBackend end

"""    MakieBackend() — dispatch tag selecting the `Makie.jl` backend."""
struct MakieBackend <: AbstractPlotBackend end

# Backends register themselves from their extension's `__init__`; the most
# recently loaded one is the default when `backend` is not given explicitly.
const _LOADED_PLOT_BACKENDS = AbstractPlotBackend[]

function _register_plot_backend(b::AbstractPlotBackend)
    b in _LOADED_PLOT_BACKENDS || push!(_LOADED_PLOT_BACKENDS, b)
    return nothing
end

"""
    default_plot_backend() -> AbstractPlotBackend

The backend used by `plot_lattice` when `backend` is not given: the most
recently loaded of `Plots` / `Makie`. Errors if neither is loaded.
"""
function default_plot_backend()
    isempty(_LOADED_PLOT_BACKENDS) && error(
        "no plotting backend loaded — run `using Plots` or `using CairoMakie` first, " *
        "or pass `backend = PlotsBackend()` / `MakieBackend()`",
    )
    return _LOADED_PLOT_BACKENDS[end]
end

_as_backend(b::AbstractPlotBackend) = b
function _as_backend(s::Symbol)
    s === :plots && return PlotsBackend()
    s === :makie && return MakieBackend()
    throw(ArgumentError("unknown plot backend :$s (use :plots or :makie)"))
end

"""
    plot_bonds!(p, lat::AbstractLattice; kwargs...)

Overlay the lattice bonds on an existing `Plots.Plot` `p` using the
per-bond wrapped displacement vector stored on each `Bond`. This
avoids the boundary-crossing "long line" artefacts that afflict naive
`position(lat, j) - position(lat, i)` plots on periodic lattices.

Concrete methods live in the Plots extension.
"""
function plot_bonds! end

"""
    plot_sites!(p, lat::AbstractLattice; kwargs...)

Overlay the lattice sites on an existing `Plots.Plot` `p`. When
`color = :auto` and the lattice has more than one geometric
sublattice, sites are grouped and coloured by sublattice id.

Concrete methods live in the Plots extension.
"""
function plot_sites! end

"""
    diffraction_pattern(ml::AbstractMomentumLattice; kwargs...)
    diffraction_pattern(lat::AbstractLattice, state, ml::AbstractMomentumLattice; kwargs...)

Visualise the Fourier-space "diffraction pattern". Two argument
shapes are supported:

1. **Pre-computed peaks.** Pass an
   [`AbstractMomentumLattice`](@ref) (typically a
   [`BraggPeakSet`](@ref) from `QuasiCrystal.jl`) — intensities are
   read directly from the momentum-lattice object.

2. **From a state on a real-space lattice.** Pass a real-space
   `lat` together with a `state` vector and a momentum-lattice
   `ml` describing the q-grid; intensities are computed from
   `S(q) = (1/N) |Σ_i s_i e^{-i q · r_i}|²` via
   [`structure_factor`](@ref). The fast paths from
   `LatticeCoreFFTWExt` / `LatticeCoreNFFTExt` are used when the
   matching extension is loaded.

Concrete methods live in the Plots extension and dispatch on
`dimension(ml)`:

- `DPhys = 1` → stem plot of `intensity(k)` vs `k`
- `DPhys = 2` → heatmap (regular meshes from `(lat, state, ml)`) or
  scatter (irregular peak sets) of `intensity(kₓ, k_y)`

Typical usage on a quasicrystal:

```julia
using Plots, LatticeCore, QuasiCrystal
qc = generate_fibonacci_projection(20)
peaks = bragg_peaks(qc; kmax = 20.0, intensity_cutoff = 1e-3)
diffraction_pattern(peaks)
```

Typical usage on a finite Bravais lattice with a snapshot:

```julia
using Plots, LatticeCore
lat = SimpleSquareLattice(8, 8, PeriodicAxis())
state = ones(Int8, num_sites(lat))
diffraction_pattern(lat, state, reciprocal_lattice(lat))
```

# Common keyword arguments

- `title` — plot title (default: "Diffraction pattern")
- `intensity::Symbol` — `:|S|²` (default) or `:S`. Selects between
  the squared modulus (intensity) and the raw structure factor.
- `marker_scale::Real` — multiplier applied to marker size in 2D
  scatter (default picks something reasonable; tune for your cutoff)
- `color` — `Plots.jl` colour value or gradient for 2D scatter / heatmap
- `log_intensity::Bool` — if `true`, plot `log10(I/I_max)` instead
  of `I`. Makes weaker peaks visible alongside Γ.
- Any other keyword is forwarded to the underlying `Plots.plot`.
"""
function diffraction_pattern end

# ---- Makie-specific viz stubs (methods live in LatticeCoreMakieExt) -
#
# `plot_lattice(lat; backend = :makie)` covers the shared lattice drawing.
# These two are Makie-only extras with no counterpart in the Plots foundation
# (a per-site field scatter and an S(k) heatmap).

"""
    makie_state(lat::AbstractLattice, state::AbstractVector; colormap=:RdBu,
                arrows=false, markersize=15, kwargs...) -> Makie.Figure

Per-site `state` (`length == num_sites(lat)`) as a colour-mapped scatter with a
colour bar. With `arrows = true` the entries are read as in-plane angles (XY
spins) and drawn as unit arrows. Concrete method in `LatticeCoreMakieExt`.
"""
function makie_state end

"""
    makie_structure_factor(lat::AbstractLattice, state::AbstractVector;
                           k_range=(-π, π), resolution=200, colormap=:viridis,
                           kwargs...) -> Makie.Figure

Heatmap of `S(k) = |Σ_j state_j e^{-i k·r_j}|² / N` over a `resolution ×
resolution` grid of `k = (kx, ky)`. Concrete method in `LatticeCoreMakieExt`.
"""
function makie_structure_factor end
