# Edge / bulk partition of a lattice, from the declared connectivity alone
# (no geometric thresholds). Generic over `AbstractLattice`. Supports
# symmetry-protected-topological analyses that need an explicit edge/bulk split.

# Under-coordinated boundary sites: coordination below the bulk coordination of
# their own sublattice (so a low-coordination *bulk* site on e.g. a Lieb/dice
# lattice is not mistaken for an edge site). A fully periodic sample has none.
function _boundary_sites(lat::AbstractLattice)
    N = num_sites(lat)
    N == 0 && return Int[]
    deg = [length(neighbors(lat, i)) for i in 1:N]
    bulk_deg = Dict{Int,Int}()
    @inbounds for i in 1:N
        s = sublattice(lat, i)
        bulk_deg[s] = haskey(bulk_deg, s) ? max(bulk_deg[s], deg[i]) : deg[i]
    end
    return [i for i in 1:N if deg[i] < bulk_deg[sublattice(lat, i)]]
end

"""
    edge_sites(lat::AbstractLattice; depth::Int=1) -> Vector{Int}

Sorted site indices of the **edge** region: sites within graph distance
`depth - 1` of an under-coordinated boundary site (coordination below the bulk
coordination of its sublattice). `depth = 1` is the boundary ring; each
increment peels one more layer inward. A fully periodic sample has an empty
edge. `depth ≥ 1`.

See also [`bulk_sites`](@ref), [`edge_bonds`](@ref).
"""
function edge_sites(lat::AbstractLattice; depth::Int=1)
    depth ≥ 1 || throw(ArgumentError("depth must be ≥ 1, got $depth"))
    N = num_sites(lat)
    N == 0 && return Int[]
    ring = _boundary_sites(lat)
    isempty(ring) && return Int[]
    dist = fill(-1, N)
    queue = Int[]
    for i in ring
        dist[i] = 0
        push!(queue, i)
    end
    head = 1
    while head ≤ length(queue)
        u = queue[head]
        head += 1
        dist[u] == depth - 1 && continue
        for v in neighbors(lat, u)
            if dist[v] == -1
                dist[v] = dist[u] + 1
                push!(queue, v)
            end
        end
    end
    return findall(d -> 0 ≤ d ≤ depth - 1, dist)
end

"""
    bulk_sites(lat::AbstractLattice; depth::Int=1) -> Vector{Int}

Sorted site indices in the **bulk** region — the complement of
[`edge_sites`](@ref)`(lat; depth)`. On a fully periodic sample this is every
site. `depth ≥ 1`.
"""
function bulk_sites(lat::AbstractLattice; depth::Int=1)
    edge = Set(edge_sites(lat; depth=depth))
    return [i for i in 1:num_sites(lat) if i ∉ edge]
end

"""
    edge_bonds(lat::AbstractLattice; depth::Int=1) -> Vector{Bond}

Bonds incident to at least one edge site (see [`edge_sites`](@ref)) — the bonds
touching the boundary region. Empty on a fully periodic sample. `depth ≥ 1`.
"""
function edge_bonds(lat::AbstractLattice; depth::Int=1)
    edge = Set(edge_sites(lat; depth=depth))
    isempty(edge) && return Bond{_spatial_dim(lat),_position_eltype(lat)}[]
    return [b for b in collect(bonds(lat)) if b.i ∈ edge || b.j ∈ edge]
end

_position_eltype(lat::AbstractLattice{D,T}) where {D,T} = T
