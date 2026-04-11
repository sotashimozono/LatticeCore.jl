"""
    AbstractLatticeElement

Abstract supertype identifying which geometric element of a lattice
a given degree of freedom lives on. The default for every
[`AbstractSiteType`](@ref) is [`VertexCenter`](@ref); other
`element_type` overrides let the interface describe bond / plaquette /
cell-centered variables such as dimers, gauge links, and flux
variables without breaking existing site-centered code.

See `dev/note/04_architecture/04_site_type/README.md` for the
two-approach strategy (line-graph versus multi-layer) that uses this
trait.
"""
abstract type AbstractLatticeElement end

"""Vertex-centered element: a site in the usual sense (default)."""
struct VertexCenter <: AbstractLatticeElement end

"""Bond-centered element: the midpoint of an edge."""
struct BondCenter <: AbstractLatticeElement end

"""Plaquette-centered element: the centre of a face."""
struct PlaquetteCenter <: AbstractLatticeElement end

"""Cell-centered element: the centre of a 3D (or higher) cell."""
struct CellCenter <: AbstractLatticeElement end
