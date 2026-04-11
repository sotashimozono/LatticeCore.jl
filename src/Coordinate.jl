"""
    AbstractCoordinate{D}

Abstract supertype for lattice coordinate representations. Concrete
subtypes describe the same lattice point in different spaces:

- [`RealSpace`](@ref) — Cartesian coordinates
- [`LatticeCoord`](@ref) — unit cell index + **geometric** sublattice
- [`HigherDimCoord`](@ref) — higher-dimensional projection (for
  cut-and-project quasicrystals)

`D` is the **physical** dimension (the dimension of the real space
the lattice lives in), which is the same for every concrete subtype
that describes a given lattice point.

Conversions between coordinate systems are performed by
[`to_real`](@ref), [`to_lattice`](@ref), and [`to_hyper`](@ref).
These functions dispatch on the concrete lattice type because the
conversion depends on the lattice basis and any higher-dimensional
projection the lattice uses.
"""
abstract type AbstractCoordinate{D} end

# ---- RealSpace -------------------------------------------------------

"""
    RealSpace{D, T}(x::SVector{D, T})
    RealSpace(x::NTuple{D, T})

Real-space Cartesian coordinate in `D` dimensions with element type
`T`.
"""
struct RealSpace{D,T<:Real} <: AbstractCoordinate{D}
    x::SVector{D,T}
end

# The `RealSpace(x::SVector{D, T})` outer constructor is auto-generated
# by Julia from the inner; we only need an extra tuple-based overload.
RealSpace(xs::NTuple{D,T}) where {D,T<:Real} = RealSpace(SVector(xs))

# ---- LatticeCoord ----------------------------------------------------

"""
    LatticeCoord{D}(cell::NTuple{D, Int}, sublattice::Int = 1)

Lattice coordinate: per-axis unit cell index plus the (1-based)
**geometric** sublattice id within the cell.

The `sublattice` field refers strictly to the *geometric* sublattice
(the honeycomb A/B positions, the Kagome A/B/C, ...), not to the
physical `AbstractSiteType` living on the site. Site types live on a
separate axis and are designed in the site-type chapter of the
architecture notes.
"""
struct LatticeCoord{D} <: AbstractCoordinate{D}
    cell::NTuple{D,Int}
    sublattice::Int

    function LatticeCoord{D}(cell::NTuple{D,Int}, sublattice::Int=1) where {D}
        return new{D}(cell, sublattice)
    end
end

function LatticeCoord(cell::NTuple{D,Int}, sublattice::Int=1) where {D}
    LatticeCoord{D}(cell, sublattice)
end

# ---- HigherDimCoord --------------------------------------------------

"""
    HigherDimCoord{DPhys, DHyper, T}(hyper::SVector{DHyper, T})

Higher-dimensional coordinate used by cut-and-project quasicrystals.
`DPhys` is the physical dimension (the target of the projection) and
`DHyper` is the host dimension (`DHyper > DPhys`). The lattice's
projection matrix maps `hyper` back into `DPhys`-dimensional real
space.
"""
struct HigherDimCoord{DPhys,DHyper,T<:Real} <: AbstractCoordinate{DPhys}
    hyper::SVector{DHyper,T}
end

function HigherDimCoord{DPhys}(hyper::SVector{DHyper,T}) where {DPhys,DHyper,T<:Real}
    return HigherDimCoord{DPhys,DHyper,T}(hyper)
end

# ---- Conversion generic functions -----------------------------------

"""
    to_real(lat::AbstractLattice, coord::AbstractCoordinate) → RealSpace

Convert `coord` into a real-space Cartesian coordinate on `lat`.
Concrete lattices should implement at least `to_real(lat, ::LatticeCoord)`
(and, for quasicrystals, `to_real(lat, ::HigherDimCoord)`).
"""
function to_real end

"""
    to_lattice(lat::AbstractLattice, coord::AbstractCoordinate) → LatticeCoord

Convert `coord` into a lattice coordinate on `lat`. Concrete lattices
should implement at least `to_lattice(lat, ::RealSpace)`.
"""
function to_lattice end

"""
    to_hyper(lat::AbstractLattice, coord::AbstractCoordinate) → HigherDimCoord

Convert `coord` into a higher-dimensional hyper coordinate. Only
meaningful for cut-and-project quasicrystals; other lattices may
leave this unimplemented.
"""
function to_hyper end

# Identity conversions — converting a coordinate to its own kind is a
# no-op regardless of the lattice. These let downstream code write a
# single `to_real(lat, coord)` call without a type check.
to_real(::Any, r::RealSpace) = r
to_lattice(::Any, lc::LatticeCoord) = lc
to_hyper(::Any, h::HigherDimCoord) = h
