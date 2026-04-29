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

Calling `plot_lattice` without loading `Plots` first raises a
`MethodError` that directs the user to the extension.
"""
function plot_lattice end

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
