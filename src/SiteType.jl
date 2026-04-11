"""
    AbstractSiteType

Abstract supertype for site types — value-level descriptors of the
physical degree of freedom living on a lattice site.

Concrete subtypes are deliberately lightweight (typically singletons
or thin parametric singletons) so they can be embedded in lattice
structs and used as a dispatch key by MC models.

# Required interface

- `state_type(st)::Type` — the Julia type used to store a state
- `random_state(rng, st)` — sample a uniformly random state

# Optional interface

- `zero_state(st)` — a canonical zero state (if defined)
- `domain(st)` — iterable over the state space (for small discrete
  site types)
- `element_type(st)::AbstractLatticeElement` — which geometric element
  the DOF lives on. Defaults to [`VertexCenter`](@ref). Override for
  bond / plaquette-centered variables such as dimer or gauge fields.

See `dev/note/04_architecture/04_site_type/README.md`.
"""
abstract type AbstractSiteType end

# ---- Required / optional generic functions ----

"""
    state_type(st::AbstractSiteType) → Type

Julia type used to store a state of site type `st`.
"""
function state_type end

"""
    random_state(rng, st::AbstractSiteType)

Sample a uniformly random state on site type `st`.
"""
function random_state end

"""
    zero_state(st::AbstractSiteType)

Canonical zero state (e.g. `0` for Ising, `0.0` for XY). Not all site
types have to define this.
"""
function zero_state end

"""
    domain(st::AbstractSiteType)

Iterable over the state space. Only meaningful for small discrete
site types (Ising, Potts); continuous types leave this unimplemented.
"""
function domain end

"""
    element_type(st::AbstractSiteType) → AbstractLatticeElement

Which geometric element the DOF lives on. Defaults to
[`VertexCenter`](@ref). Override for bond / plaquette / cell-centered
site types.
"""
element_type(::AbstractSiteType) = VertexCenter()

# ---- Built-in site types ---------------------------------------------

"""
    IsingSite{T <: Integer}()

Ising spin: state is `±1` stored as `T` (default `Int8`).
"""
struct IsingSite{T<:Integer} <: AbstractSiteType end
IsingSite() = IsingSite{Int8}()

state_type(::IsingSite{T}) where {T} = T
random_state(rng, ::IsingSite{T}) where {T} = rand(rng, (T(-1), T(1)))
zero_state(::IsingSite{T}) where {T} = zero(T)
domain(::IsingSite{T}) where {T} = (T(-1), T(1))

"""
    PottsSite{Q, T <: Integer}()

`Q`-state Potts spin: state is an integer in `1..Q` stored as `T`.
"""
struct PottsSite{Q,T<:Integer} <: AbstractSiteType end
PottsSite(q::Integer) = PottsSite{Int(q),Int8}()

state_type(::PottsSite{Q,T}) where {Q,T} = T
random_state(rng, ::PottsSite{Q,T}) where {Q,T} = T(rand(rng, 1:Q))
domain(::PottsSite{Q,T}) where {Q,T} = ntuple(i -> T(i), Q)

"""
    XYSite{T <: AbstractFloat}()

Planar rotor: state is an angle `θ ∈ [0, 2π)` stored as `T`
(default `Float64`).
"""
struct XYSite{T<:AbstractFloat} <: AbstractSiteType end
XYSite() = XYSite{Float64}()

state_type(::XYSite{T}) where {T} = T
random_state(rng, ::XYSite{T}) where {T} = T(2π) * rand(rng, T)
zero_state(::XYSite{T}) where {T} = zero(T)

"""
    HeisenbergSite{T <: AbstractFloat}()

Classical Heisenberg spin: state is a unit vector in ℝ³ stored as
`SVector{3, T}` (default `Float64`).
"""
struct HeisenbergSite{T<:AbstractFloat} <: AbstractSiteType end
HeisenbergSite() = HeisenbergSite{Float64}()

state_type(::HeisenbergSite{T}) where {T} = SVector{3,T}
function random_state(rng, ::HeisenbergSite{T}) where {T}
    # Marsaglia-style uniform point on S².
    z = T(2) * rand(rng, T) - one(T)
    φ = T(2π) * rand(rng, T)
    sxy = sqrt(one(T) - z * z)
    return SVector{3,T}(sxy * cos(φ), sxy * sin(φ), z)
end

"""
    EmptySite()

Vacancy / empty site: stores `nothing`. Useful for diluted models
and as a placeholder where a site has no degrees of freedom.
"""
struct EmptySite <: AbstractSiteType end

state_type(::EmptySite) = Nothing
random_state(::Any, ::EmptySite) = nothing
zero_state(::EmptySite) = nothing
