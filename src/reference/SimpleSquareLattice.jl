"""
    SimpleSquareLattice{T, B <: AbstractBoundaryCondition}(Lx, Ly, boundary)

A 2D square lattice of `Lx × Ly` sites with unit spacing. The type
parameter `B` selects the boundary condition (`PBC` or `OBC`), applied
uniformly to both axes.

This is LatticeCore's 2D reference implementation. It exists so that
the core interface can be exercised end-to-end without depending on
`Lattice2D.jl`. The production-grade 2D lattice, with the full
topology catalogue (triangular, honeycomb, kagome, ...), sublattice
site types, and reciprocal lattice machinery, lives in `Lattice2D.jl`.

# Site indexing

Sites are laid out in row-major order:

    site_index(x, y) = (y - 1) * Lx + x     # 1 <= x <= Lx, 1 <= y <= Ly

so sites 1..Lx form the bottom row, Lx+1..2Lx the next, and so on.

# Example
```julia
julia> lat = SimpleSquareLattice(3, 3, PBC());

julia> num_sites(lat)
9

julia> position(lat, 5)   # middle of a 3x3 lattice
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 2.0
 2.0
```
"""
struct SimpleSquareLattice{T<:AbstractFloat,B<:AbstractBoundaryCondition} <:
       AbstractLattice{2,T}
    Lx::Int
    Ly::Int
    boundary::B
end

# Convenience constructor: default to Float64 positions and PBC.
function SimpleSquareLattice(
    Lx::Int, Ly::Int, boundary::B=PBC()
) where {B<:AbstractBoundaryCondition}
    return SimpleSquareLattice{Float64,B}(Lx, Ly, boundary)
end

# ---- Private row-major coordinate helpers ----

@inline function _site_to_xy(l::SimpleSquareLattice, i::Int)
    x = mod1(i, l.Lx)
    y = ((i - 1) ÷ l.Lx) + 1
    return (x, y)
end

@inline function _xy_to_site(l::SimpleSquareLattice, x::Int, y::Int)
    return (y - 1) * l.Lx + x
end

# ---- Required interface ----

num_sites(l::SimpleSquareLattice) = l.Lx * l.Ly

function position(l::SimpleSquareLattice{T}, i::Int) where {T}
    x, y = _site_to_xy(l, i)
    return SVector{2,T}(T(x), T(y))
end

function neighbors(l::SimpleSquareLattice{T,PBC}, i::Int) where {T}
    x, y = _site_to_xy(l, i)
    ns = [
        _xy_to_site(l, mod1(x + 1, l.Lx), y),
        _xy_to_site(l, mod1(x - 1, l.Lx), y),
        _xy_to_site(l, x, mod1(y + 1, l.Ly)),
        _xy_to_site(l, x, mod1(y - 1, l.Ly)),
    ]
    # Collapse degenerate Lx == 2 / Ly == 2 cases that would otherwise
    # return duplicate neighbours.
    return unique!(ns)
end

function neighbors(l::SimpleSquareLattice{T,OBC}, i::Int) where {T}
    x, y = _site_to_xy(l, i)
    ns = Int[]
    if x + 1 <= l.Lx
        push!(ns, _xy_to_site(l, x + 1, y))
    end
    if x - 1 >= 1
        push!(ns, _xy_to_site(l, x - 1, y))
    end
    if y + 1 <= l.Ly
        push!(ns, _xy_to_site(l, x, y + 1))
    end
    if y - 1 >= 1
        push!(ns, _xy_to_site(l, x, y - 1))
    end
    return ns
end

boundary(l::SimpleSquareLattice) = l.boundary

size_trait(l::SimpleSquareLattice) = FiniteSize((l.Lx, l.Ly))

# ---- Trait overrides ----

topology(::SimpleSquareLattice) = TopologyTrait{:square}()

periodicity(::SimpleSquareLattice{T,PBC}) where {T} = Periodic()
periodicity(::SimpleSquareLattice{T,OBC}) where {T} = Aperiodic()

# PBC square lattice is bipartite iff both axis lengths are even.
function is_bipartite(l::SimpleSquareLattice{T,PBC}) where {T}
    return iseven(l.Lx) && iseven(l.Ly)
end
is_bipartite(::SimpleSquareLattice{T,OBC}) where {T} = true

reciprocal_support(::SimpleSquareLattice{T,PBC}) where {T} = HasReciprocal()
reciprocal_support(::SimpleSquareLattice{T,OBC}) where {T} = NoReciprocal()
