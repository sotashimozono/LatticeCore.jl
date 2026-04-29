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

Sine-square deformation modifier. `L` is the characteristic scale
used by the envelope `sin²(π · x/L)`. The concrete weight implementation
will land with the MC layer; for now this serves as a placeholder
carrying the scale parameter.

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
# SSD weight is intentionally left unimplemented for now; it will land
# with the MC layer once we have a `center(lat)` accessor.
