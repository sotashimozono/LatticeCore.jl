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
"""
function bond_center(lat::AbstractLattice{D,T}, bond::Bond{D,T}) where {D,T}
    return (position(lat, bond.i) + position(lat, bond.j)) / 2
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
