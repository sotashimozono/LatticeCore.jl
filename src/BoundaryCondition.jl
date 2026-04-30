"""
    AbstractBoundaryCondition

Abstract supertype for lattice boundary conditions.

The canonical concrete type is [`LatticeBoundary`](@ref), which
composes:

- a tuple of per-axis boundary conditions (subtypes of
  [`AbstractAxisBC`](@ref): [`PeriodicAxis`](@ref),
  [`OpenAxis`](@ref), [`TwistedAxis`](@ref))
- a non-topological modifier (subtype of
  [`AbstractBoundaryModifier`](@ref): [`NoModifier`](@ref),
  [`SSD`](@ref))

The split is deliberate: an *axis BC* decides whether a candidate
bond exists (and whether it carries a twist phase), while a *modifier*
only reweights existing bonds (e.g. sine-square deformation). Multiple
concerns are composed, not shoehorned into a single taxonomy.

See `dev/note/04_architecture/03_boundary_and_coordinates/README.md`
for the design rationale.
"""
abstract type AbstractBoundaryCondition end

# ---- Axis-level BC ---------------------------------------------------

"""
    AbstractAxisBC

Abstract supertype for per-axis boundary conditions. Concrete subtypes
decide how a raw candidate cell index is wrapped into the lattice's
index range along a single axis.

Subtypes:
- [`PeriodicAxis`](@ref) — wraps via `mod1`
- [`OpenAxis`](@ref) — rejects out-of-range indices
- [`TwistedAxis`](@ref) — wraps like `PeriodicAxis` but attaches a
  phase factor to bonds that cross the boundary
"""
abstract type AbstractAxisBC end

"""Periodic axis BC: wrap-around via `mod1`."""
struct PeriodicAxis <: AbstractAxisBC end

"""Open axis BC: bonds crossing the boundary are dropped."""
struct OpenAxis <: AbstractAxisBC end

"""
    TwistedAxis(phase::Real)

Periodic axis BC with a phase angle (radians) attached to
boundary-crossing bonds. The connectivity is identical to
[`PeriodicAxis`](@ref); only the phase returned by [`axis_phase`](@ref)
differs.
"""
struct TwistedAxis{T<:Real} <: AbstractAxisBC
    phase::T
end

# ---- Modifier (non-topological) --------------------------------------

"""
    AbstractBoundaryModifier

Abstract supertype for boundary modifiers. A modifier does not change
which bonds exist — it only reweights existing bonds through
[`bond_weight`](@ref). Examples: [`NoModifier`](@ref),
[`SSD`](@ref).

# Extending

Downstream packages may define their own modifier types by
subtyping `AbstractBoundaryModifier` and overloading

    bond_weight(::MyModifier, lat::AbstractLattice, i::Int, j::Int)::Float64

The two-argument default `bond_weight(modifier, lat)` may also be
overloaded if a modifier wants to broadcast a single scalar over
the whole lattice. `AbstractBoundaryModifier` and the concrete
[`NoModifier`](@ref) / [`SSD`](@ref) types are part of the public
exported API of `LatticeCore` precisely so downstream code can
extend them.
"""
abstract type AbstractBoundaryModifier end

"""Identity modifier: every bond has weight `1.0`."""
struct NoModifier <: AbstractBoundaryModifier end

"""
    SSD(L)

Sine-square deformation modifier. `L` is a characteristic scale
carried by the modifier; the canonical bond weight evaluation
([`bond_weight(::SSD, lat, i, j)`](@ref bond_weight)) reads the
per-axis lengths of the lattice from [`size_trait`](@ref) instead, so
`L` is informational under finite, fully-specified geometries. It is
retained for downstream packages that may want to override the
default with an alternative scale (e.g. infinite-system extrapolations).

# Canonical envelope

For a `D`-dimensional finite lattice with per-axis cell counts
`(L_1, …, L_D)`, the SSD envelope evaluated at a lattice cell with
1-based per-axis coordinate `cx_d ∈ 1:L_d` is

```math
f(\\mathbf{r}) = \\prod_{d=1}^{D} \\sin^{2}\\!\\left(
    \\pi \\, \\frac{c_{x,d} - 1/2}{L_d}
\\right),
```

equivalent to the standard ``\\sin^2(\\pi (i + 1/2) / L)`` on
0-indexed sites. The bond weight between sites `i, j` is the
arithmetic mean of the two endpoint envelopes,
``w(i, j) = (f(\\mathbf{r}_i) + f(\\mathbf{r}_j)) / 2``.

`SSD` is exported so downstream MC / TN packages can construct
boundaries with SSD weighting and dispatch on the type. Custom
deformations should subtype [`AbstractBoundaryModifier`](@ref)
rather than re-using `SSD`.
"""
struct SSD{T<:Real} <: AbstractBoundaryModifier
    L::T
end

# ---- Composite boundary ----------------------------------------------

"""
    LatticeBoundary{N, A, M}(axes::NTuple{N, <:AbstractAxisBC}, modifier)

Composite boundary condition for an `N`-dimensional lattice. Stores
one [`AbstractAxisBC`](@ref) per axis plus a single
[`AbstractBoundaryModifier`](@ref). Supports mixed-axis boundary
conditions natively — a cylinder is

```julia
LatticeBoundary((PeriodicAxis(), OpenAxis()))
```

The zero-modifier form is the default.
"""
struct LatticeBoundary{N,A<:Tuple{Vararg{AbstractAxisBC}},M<:AbstractBoundaryModifier} <:
       AbstractBoundaryCondition
    axes::A
    modifier::M

    function LatticeBoundary(
        axes::A, modifier::M
    ) where {A<:Tuple{Vararg{AbstractAxisBC}},M<:AbstractBoundaryModifier}
        return new{length(axes),A,M}(axes, modifier)
    end
end

"""
    LatticeBoundary(axes::NTuple{N, <:AbstractAxisBC})

Convenience constructor: composite boundary with [`NoModifier`](@ref).
"""
LatticeBoundary(axes::Tuple{Vararg{AbstractAxisBC}}) = LatticeBoundary(axes, NoModifier())

# ---- Hook functions --------------------------------------------------

"""
    apply_axis_bc(axis_bc::AbstractAxisBC, idx::Int, L::Int) → (wrapped, is_valid)

Apply the axis-level BC to a raw cell index `idx ∈ ℤ` on an axis of
length `L`. Returns a tuple `(wrapped::Int, is_valid::Bool)`:

- `PeriodicAxis`, `TwistedAxis`: `(mod1(idx, L), true)`
- `OpenAxis`: `(idx, 1 <= idx <= L)`

The twist *phase* is intentionally reported separately by
[`axis_phase`](@ref): MC code paths that do not need phases (classical
Ising / XY / Heisenberg) are not forced to traffic in complex numbers.
"""
apply_axis_bc(::PeriodicAxis, idx::Int, L::Int) = (mod1(idx, L), true)
apply_axis_bc(::OpenAxis, idx::Int, L::Int) = (idx, 1 <= idx <= L)
apply_axis_bc(::TwistedAxis, idx::Int, L::Int) = (mod1(idx, L), true)

"""
    axis_phase(axis_bc::AbstractAxisBC, idx::Int, L::Int) → ComplexF64

Phase factor attached to a bond that steps from a valid cell index
to `idx` (possibly out of `1:L`) along a single axis.

- `PeriodicAxis`, `OpenAxis`: always `1 + 0im`.
- `TwistedAxis(θ)`: `cis(+θ)` if `idx > L`, `cis(-θ)` if `idx < 1`, otherwise `1 + 0im`.

Classical MC paths may ignore this; quantum / flux-carrying code
should multiply bond contributions by the phase.
"""
axis_phase(::PeriodicAxis, ::Int, ::Int) = complex(1.0, 0.0)
axis_phase(::OpenAxis, ::Int, ::Int) = complex(1.0, 0.0)
function axis_phase(ax::TwistedAxis, idx::Int, L::Int)
    if idx > L
        return cis(float(ax.phase))
    elseif idx < 1
        return cis(-float(ax.phase))
    else
        return complex(1.0, 0.0)
    end
end

"""
    bond_weight(modifier::AbstractBoundaryModifier, lat, i::Int, j::Int) → Float64

Multiplicative weight applied to the bond `(i, j)` by a boundary
modifier. Default implementation for [`NoModifier`](@ref) returns
`1.0`. SSD and other non-trivial modifiers will override this.
"""
bond_weight(::NoModifier, lat, ::Int, ::Int) = 1.0

"""
    _ssd_axis_envelope(cx::Real, L::Integer) → Float64

Single-axis SSD envelope `sin²(π (cx - 1/2) / L)` for a 1-based
real-valued per-axis cell coordinate `cx ∈ [1, L]`. The shift by
`-1/2` converts the convention to the canonical 0-indexed form
`sin²(π (i + 1/2) / L)` with `i = cx - 1`.

Internal helper for [`bond_weight(::SSD, lat, i, j)`](@ref bond_weight).
"""
@inline function _ssd_axis_envelope(cx::Real, L::Integer)
    return sin(pi * (float(cx) - 0.5) / L)^2
end

"""
    _ssd_site_envelope(lat::AbstractLattice, i::Int) → Float64

Multi-axis SSD envelope evaluated at site `i`: the product of
[`_ssd_axis_envelope`](@ref) over the per-axis cell coordinates of
`i`. Uses [`size_trait`](@ref) for per-axis lengths and
[`to_lattice`](@ref) ∘ [`position`](@ref) for the cell coordinate
of the site.

The conversion `to_lattice(lat, RealSpace(position(lat, i)))` round-trips
exactly for the reference lattices (which use unit spacing) and relies
on the lattice's own `to_lattice(::RealSpace)` for any custom basis.
"""
@inline function _ssd_site_envelope(lat::AbstractLattice, i::Int)
    st = size_trait(lat)
    st isa FiniteSize ||
        throw(ArgumentError("SSD weight requires a FiniteSize lattice; got $(typeof(st))"))
    dims = st.dims
    coord = to_lattice(lat, RealSpace(position(lat, i)))
    f = 1.0
    @inbounds for d in eachindex(dims)
        f *= _ssd_axis_envelope(coord.cell[d], dims[d])
    end
    return f
end

"""
    bond_weight(::SSD, lat::AbstractLattice, i::Int, j::Int) → Float64

Sine-square deformation envelope for the bond `(i, j)`. Each axis
contributes a `sin²(π (c_d - 1/2) / L_d)` factor, where `c_d` is the
1-based per-axis cell coordinate read from
`to_lattice(lat, RealSpace(position(lat, i)))` and `L_d` is the
matching component of `size_trait(lat).dims`. The bond weight is the
arithmetic mean of the two endpoint envelopes,

```math
w(i, j) = \\tfrac{1}{2}\\bigl(f(\\mathbf{r}_i) + f(\\mathbf{r}_j)\\bigr).
```

Requires a [`FiniteSize`](@ref) lattice; an
[`InfiniteSize`](@ref) / [`QuasiInfiniteSize`](@ref) lattice raises
`ArgumentError` because the envelope has no canonical scale without
finite per-axis lengths. Downstream packages may overload the method
on a more specific lattice type to supply alternative scales.

# Boundary-condition assumptions

SSD is most commonly paired with open boundaries (the envelope
suppresses surface excitations of an OBC sample); the implementation
itself is BC-agnostic and will weight any bond regardless of the axis
BCs in [`LatticeBoundary`](@ref).
"""
function bond_weight(modifier::SSD, lat::AbstractLattice, i::Int, j::Int)
    fi = _ssd_site_envelope(lat, i)
    fj = _ssd_site_envelope(lat, j)
    return 0.5 * (fi + fj)
end
