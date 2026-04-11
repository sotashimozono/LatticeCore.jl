"""
    AbstractSiteLayout

Abstract supertype for strategies that store the per-site
[`AbstractSiteType`](@ref) of a lattice. Three concrete layouts are
provided, each picking a different memory / flexibility trade-off:

- [`UniformLayout`](@ref) — every site shares the same site type
  (stored once; zero per-site memory overhead)
- [`SublatticeLayout`](@ref) — each geometric sublattice has its own
  site type (`Vector{Int}` of sublattice ids)
- [`ExplicitLayout`](@ref) — one `AbstractSiteType` per site (for
  disordered / quenched-random configurations)

All layouts satisfy `site_type(layout, i)::AbstractSiteType`.
"""
abstract type AbstractSiteLayout end

"""
    site_type(layout::AbstractSiteLayout, i::Int) → AbstractSiteType

Return the site type stored at site index `i`.
"""
function site_type end

# ---- UniformLayout ---------------------------------------------------

"""
    UniformLayout(st::AbstractSiteType)

Every site has the same site type. The type parameter carries the
concrete site-type singleton so that downstream MC code can dispatch
on it at compile time (constant folding in the hot loop).
"""
struct UniformLayout{S<:AbstractSiteType} <: AbstractSiteLayout
    site_type::S
end

site_type(l::UniformLayout, ::Int) = l.site_type

# ---- SublatticeLayout ------------------------------------------------

"""
    SublatticeLayout(by_sublattice::NTuple{N, <:AbstractSiteType}, sublattice_of::Vector{Int})

Per-sublattice site type. `by_sublattice` is a tuple of site type
instances, one per geometric sublattice; `sublattice_of[i]` tells
which sublattice site `i` belongs to (1-based).

This is the natural layout for mixed-spin models (e.g. Ising on the
A sublattice, XY on B).
"""
struct SublatticeLayout{N,Tup<:NTuple{N,AbstractSiteType}} <: AbstractSiteLayout
    by_sublattice::Tup
    sublattice_of::Vector{Int}

    function SublatticeLayout(
        by_sublattice::Tup, sublattice_of::Vector{Int}
    ) where {N,Tup<:NTuple{N,AbstractSiteType}}
        return new{N,Tup}(by_sublattice, sublattice_of)
    end
end

site_type(l::SublatticeLayout, i::Int) = l.by_sublattice[l.sublattice_of[i]]

# ---- ExplicitLayout --------------------------------------------------

"""
    ExplicitLayout(types::Vector{<:AbstractSiteType})

One site type per site. Use for quenched disorder or any
configuration where the per-site type cannot be factored through a
sublattice assignment.
"""
struct ExplicitLayout{S<:AbstractSiteType} <: AbstractSiteLayout
    types::Vector{S}
end

site_type(l::ExplicitLayout, i::Int) = l.types[i]
