module LatticeCore

"""
    AbstractLattice{D}
abstract type for lattices in D dimensions.
"""
abstract type AbstractLattice{D} end
export AbstractLattice

"""
    FilteredLattice{D, L<:AbstractLattice{D}, F<:Function}
struct for a filtered lattice, which is a sub-lattice defined by a filtering function.
- `parent::L`: the parent lattice from which this filtered lattice is derived.
- `filter_func::F`: a function that takes a coordinate vector and returns a Bool indicating whether the site is included in the filtered lattice.
- `active_indices::Vector{Int}`: indices of the sites in the parent lattice that are included in the filtered lattice.
"""
struct FilteredLattice{D, L<:AbstractLattice{D}, F<:Function} <: AbstractLattice{D}
    parent::L
    filter_func::F # r -> norm(r) < R などの条件
    active_indices::Vector{Int}
end
"""
    AbstractLatticeConnection
abstract type for lattice connections (edges, bonds).
"""
abstract type AbstractLatticeConnection end
export AbstractLatticeConnection

"""
    Bond
struct for a bond (edge) in the lattice.
- `src::Int`: start site index
- `dst::Int`: destination site index
- `type::Int`: type of the bond for categorization
- `vector::Vector{Float64}`: expresses the bond vector from src to dst
"""
struct Bond <: AbstractLatticeConnection
    src::Int
    dst::Int
    type::Int
    vector::Vector{Float64}
end
export Bond

"""
    AbstractBoundaryCondition
Abstract type for boundary conditions.
"""
abstract type AbstractBoundaryCondition end
export AbstractBoundaryCondition

"""
    AbstractIndexing
Abstract type for indexing schemes.
"""
abstract type AbstractIndexing end
export AbstractIndexing

"""
    AbstractConstructionMethod
Abstract type for lattice construction methods.
"""
abstract type AbstractConstructionMethod end
export AbstractConstructionMethod

end
