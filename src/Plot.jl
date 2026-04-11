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

Visualise the Fourier-space "diffraction pattern" of a momentum
lattice. Concrete methods live in the Plots extension and dispatch
on the physical dimension of the momentum lattice:

- `DPhys = 1` → stem plot of `intensity(k)` vs `k`
- `DPhys = 2` → scatter of peak positions in `(kₓ, k_y)` with
  marker size (and colour) proportional to intensity

Typical usage on a quasicrystal:

```julia
using Plots, LatticeCore, QuasiCrystal
qc = generate_fibonacci_projection(20)
peaks = bragg_peaks(qc; kmax = 20.0, intensity_cutoff = 1e-3)
diffraction_pattern(peaks)
```

# Common keyword arguments

- `title` — plot title (default: "Diffraction pattern")
- `marker_scale::Real` — multiplier applied to marker size in 2D
  (default picks something reasonable; tune for your cutoff)
- `color` — `Plots.jl` colour value or gradient for 2D scatter
- `log_intensity::Bool` — if `true`, plot `log10(I/I_max)` instead
  of `I`. Makes weaker peaks visible alongside Γ.
- Any other keyword is forwarded to the underlying `Plots.plot`.
"""
function diffraction_pattern end
