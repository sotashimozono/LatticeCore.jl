"""
    LineLattice{T, B <: AbstractBoundaryCondition}(N, boundary)

A 1D linear chain of `N` sites with unit spacing. The type parameter
`B` selects the boundary condition: `PBC` (periodic) or `OBC` (open).

This is LatticeCore's simplest reference implementation; paired with
`SimpleSquareLattice` it is also the canonical mock used by downstream
Monte Carlo unit tests.

# Example
```julia
julia> lat = LineLattice(5, PBC());

julia> num_sites(lat)
5

julia> neighbors(lat, 1)
2-element Vector{Int64}:
 5
 2
```
"""
struct LineLattice{T<:AbstractFloat,B<:AbstractBoundaryCondition} <: AbstractLattice{1,T}
    N::Int
    boundary::B
end

# Convenience constructor: default to Float64 positions.
function LineLattice(N::Int, boundary::B=PBC()) where {B<:AbstractBoundaryCondition}
    return LineLattice{Float64,B}(N, boundary)
end

# ---- Required interface ----

num_sites(l::LineLattice) = l.N

position(l::LineLattice{T}, i::Int) where {T} = SVector{1,T}(T(i))

function neighbors(l::LineLattice{T,PBC}, i::Int) where {T}
    # `unique!` collapses the degenerate N == 2 case (both neighbours
    # point to the other site) to a single entry.
    return unique!([mod1(i - 1, l.N), mod1(i + 1, l.N)])
end

function neighbors(l::LineLattice{T,OBC}, i::Int) where {T}
    return filter(j -> 1 <= j <= l.N, [i - 1, i + 1])
end

boundary(l::LineLattice) = l.boundary

size_trait(l::LineLattice) = FiniteSize((l.N,))

# ---- Trait overrides ----

topology(::LineLattice) = TopologyTrait{:line}()

periodicity(::LineLattice{T,PBC}) where {T} = Periodic()
periodicity(::LineLattice{T,OBC}) where {T} = Aperiodic()

# PBC chain is bipartite iff the cycle length N is even.
is_bipartite(l::LineLattice{T,PBC}) where {T} = iseven(l.N)
is_bipartite(::LineLattice{T,OBC}) where {T} = true

reciprocal_support(::LineLattice{T,PBC}) where {T} = HasReciprocal()
reciprocal_support(::LineLattice{T,OBC}) where {T} = NoReciprocal()
