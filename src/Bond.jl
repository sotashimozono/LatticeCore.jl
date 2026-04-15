"""
    Bond{D, T}

A bond (edge) connecting two sites of a lattice.

# Fields
- `i::Int` — source site index
- `j::Int` — target site index
- `vector::SVector{D, T}` — displacement `position(j) - position(i)`,
  wrapped by the lattice's boundary condition if applicable
- `type::Symbol` — bond type tag (e.g. `:nearest`, `:next_nearest`,
  `:dimer_strong`). Used as a dispatch key by anisotropic models
  (see the 07 MC layer design note).
"""
struct Bond{D,T}
    i::Int
    j::Int
    vector::SVector{D,T}
    type::Symbol
end

"""
    bond_center(lat::AbstractLattice, bond::Bond) → SVector{D, T}

Geometric center (midpoint) of the bond in real space. Useful for
bond-centered observables and for position-dependent bond modifiers
such as sine-square deformation.

Uses the bond's stored displacement `bond.vector` rather than the
literal positions of `bond.i` and `bond.j`. This matters under
periodic boundary conditions: the wrapped target's `position(j)` is
on the opposite side of the sample, so `(position(i) + position(j))/2`
would point to the sample interior instead of just outside the
boundary. The bond carries the unwrapped displacement, so
`position(i) + bond.vector / 2` gives the true geometric midpoint.
"""
function bond_center(lat::AbstractLattice{D,T}, bond::Bond{D,T}) where {D,T}
    return position(lat, bond.i) + bond.vector / 2
end

"""
    bonds(lat::AbstractLattice)

Iterator of all bonds in the lattice. The default implementation
builds `Bond` objects on the fly from `neighbors(lat, i)` using the
`:nearest` tag. Concrete lattices may override this for efficiency or
to attach non-default bond types (e.g. dimer-strong vs dimer-weak).
"""
function bonds(lat::AbstractLattice{D,T}) where {D,T}
    return (
        Bond{D,T}(i, j, position(lat, j) - position(lat, i), :nearest) for
        i in 1:num_sites(lat) for j in neighbors(lat, i) if j > i
    )
end

"""
    neighbor_bonds(lat::AbstractLattice, i::Int)

Iterator of bonds incident to site `i`. The default implementation
builds `Bond` objects from `neighbors(lat, i)` using the `:nearest`
tag. This is the canonical entry point the 07 MC layer uses to walk
interactions involving a given site.
"""
function neighbor_bonds(lat::AbstractLattice{D,T}, i::Int) where {D,T}
    return (
        Bond{D,T}(i, j, position(lat, j) - position(lat, i), :nearest) for
        j in neighbors(lat, i)
    )
end

# ---- Generic element-center accessors --------------------------------
#
# These methods present a uniform API for "ask the lattice about its
# elements of centring X" without having to remember whether the
# concept is sites, bonds, plaquettes, or cells. Concrete lattices may
# specialise any of them for efficiency or to support `PlaquetteCenter`
# / `CellCenter`. The function declarations live in `LatticeElement.jl`;
# the methods need both `AbstractLattice` and `Bond` so they live here.

# num_elements ---------------------------------------------------------

num_elements(lat::AbstractLattice, ::VertexCenter) = num_sites(lat)
function num_elements(lat::AbstractLattice, ::BondCenter)
    return count(_ -> true, bonds(lat))
end

# elements -------------------------------------------------------------

elements(lat::AbstractLattice, ::VertexCenter) = 1:num_sites(lat)
elements(lat::AbstractLattice, ::BondCenter) = bonds(lat)

# element_position -----------------------------------------------------

element_position(lat::AbstractLattice, ::VertexCenter, i::Int) = position(lat, i)
function element_position(lat::AbstractLattice, ::BondCenter, i::Int)
    bs = collect(bonds(lat))
    return bond_center(lat, bs[i])
end

# element_positions ----------------------------------------------------

function element_positions(lat::AbstractLattice, e::AbstractLatticeElement)
    return (element_position(lat, e, i) for i in 1:num_elements(lat, e))
end
