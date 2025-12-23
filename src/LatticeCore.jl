module LatticeCore

"""
    AbstractLattice{D}
abstract type for lattices in D dimensions.
"""
abstract type AbstractLattice{D} end
export AbstractLattice

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
    Connection
Connection rules within or between unit cells.
- `src_sub`: sublattice index of the start point (1, 2, ...)
- `dst_sub`: sublattice index of the end point
- `dx`, `dy`: relative cell position of the end point (0,0 means within the same unit cell)
- `type`: type of the connection
"""
struct Connection <: AbstractLatticeConnection
    src_sub::Int
    dst_sub::Int
    dx::Int
    dy::Int
    type::Int
end
export Connection

"""
    AbstractTopology
Abstract type for lattice topologies.
"""
abstract type AbstractTopology{D} <: AbstractLattice{D} end

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

end
