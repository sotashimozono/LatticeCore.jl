"""
    Graph operations for `AbstractLattice`.

This file provides a small, dependency-light graph layer on top of the
`AbstractLattice` interface:

- [`adjacency_matrix`](@ref) — N × N adjacency as a `SparseMatrixCSC{Bool, Int}`
  (or a dense `Matrix{Bool}` when `sparse=false`),
- [`shortest_path`](@ref) — unweighted hop-count shortest path via BFS,
  returning both the distance and a reconstructed path,
- [`connected_components`](@ref) — connected components of the lattice
  graph as a `Vector{Vector{Int}}`.

Adjacency is built from the canonical `BondCenter` element iterator
(`elements(lat, BondCenter())`) so that whatever bonds a lattice
publishes are exactly what these graph operations see. This keeps the
API generic over any concrete `AbstractLattice` and consistent with
the rest of the element-centric interface.

Only `SparseArrays` (stdlib) is used; `Graphs.jl` is intentionally not
a dependency.
"""

using SparseArrays

"""
    adjacency_matrix(lat::AbstractLattice; sparse::Bool=true)

Return the N × N undirected adjacency matrix of the lattice graph,
where `N = num_sites(lat)`. The matrix is symmetric; `A[i, j] == true`
iff there is a bond between sites `i` and `j`.

Bonds are taken from `elements(lat, BondCenter())`, so concrete lattices
that override `bonds(lat)` automatically get a consistent adjacency
matrix without further work.

# Keyword arguments
- `sparse=true` (default): return a `SparseMatrixCSC{Bool, Int}`. This
  is the recommended form for large lattices.
- `sparse=false`: return a dense `Matrix{Bool}` of the same content.
  Useful for tiny test lattices and notebooks.

Self-loops (a bond with `i == j`) are skipped. Duplicate bonds in the
input are folded into a single entry.

# Example
```julia
lat = SimpleSquareLattice(3, 3, OpenAxis())
A = adjacency_matrix(lat)        # 9×9 SparseMatrixCSC{Bool, Int}
A == transpose(A)                # true
```
"""
function adjacency_matrix(lat::AbstractLattice; sparse::Bool=true)
    is_finite(lat) || throw(ArgumentError("adjacency_matrix requires a finite lattice"))
    N = num_sites(lat)
    Is = Int[]
    Js = Int[]
    seen = Set{Tuple{Int,Int}}()
    for b in elements(lat, BondCenter())
        i, j = b.i, b.j
        i == j && continue
        key = i < j ? (i, j) : (j, i)
        key in seen && continue
        push!(seen, key)
        push!(Is, i)
        push!(Js, j)
        push!(Is, j)
        push!(Js, i)
    end
    Vs = trues(length(Is))
    return _adjacency_matrix_assemble(Is, Js, Vs, N, sparse)
end

# Internal: pick sparse vs dense storage. `combine = |` folds any
# duplicated coordinate (shouldn't happen because of the `seen` guard,
# but cheap insurance) into a single `true`. The keyword argument
# shadowing happens because `adjacency_matrix(...; sparse=true)`
# binds `sparse` locally to a `Bool`, which would otherwise prevent
# us from calling `SparseArrays.sparse(...)`. We pass the flag in
# under a non-colliding name and qualify the constructor explicitly.
function _adjacency_matrix_assemble(
    Is::Vector{Int}, Js::Vector{Int}, Vs::AbstractVector{Bool}, N::Int, use_sparse::Bool
)
    if use_sparse
        return SparseArrays.sparse(Is, Js, Vs, N, N, |)
    else
        # Return a true `Matrix{Bool}` (not a `BitMatrix`) so callers
        # that dispatch on dense representations get the expected
        # element-type-preserving array.
        M = zeros(Bool, N, N)
        @inbounds for k in eachindex(Is)
            M[Is[k], Js[k]] = true
        end
        return M
    end
end

"""
    shortest_path(lat::AbstractLattice, src::Int, dst::Int)
        -> (dist::Int, path::Vector{Int})

Unweighted shortest path on the lattice graph from `src` to `dst`,
computed via breadth-first search using `neighbors(lat, ·)`.

Returns a tuple `(dist, path)` where:
- `dist` is the hop count (number of bonds) between `src` and `dst`.
  `0` if `src == dst`.
- `path` is the reconstructed sequence of site indices, starting at
  `src` and ending at `dst`. `[src]` if `src == dst`.

If `src` and `dst` are in different connected components,
`(typemax(Int), Int[])` is returned. Edges are treated as having unit
weight; for weighted shortest paths, downstream packages may add a
Dijkstra/A* layer on top.

# Example
```julia
lat = SimpleSquareLattice(3, 3, OpenAxis())
d, p = shortest_path(lat, 1, 9)
# d == 4, p == [1, 2, 3, 6, 9]   (one of several minimum paths)
```
"""
function shortest_path(lat::AbstractLattice, src::Int, dst::Int)
    is_finite(lat) || throw(ArgumentError("shortest_path requires a finite lattice"))
    N = num_sites(lat)
    (1 <= src <= N) || throw(BoundsError(lat, src))
    (1 <= dst <= N) || throw(BoundsError(lat, dst))

    if src == dst
        return (0, [src])
    end

    parent = fill(0, N)
    visited = falses(N)
    visited[src] = true
    queue = Int[src]
    head = 1
    found = false
    while head <= length(queue)
        u = queue[head]
        head += 1
        for v in neighbors(lat, u)
            if !visited[v]
                visited[v] = true
                parent[v] = u
                if v == dst
                    found = true
                    break
                end
                push!(queue, v)
            end
        end
        found && break
    end

    if !found
        return (typemax(Int), Int[])
    end

    # Reconstruct path src -> ... -> dst by walking parents backwards.
    path = Int[dst]
    cur = dst
    while cur != src
        cur = parent[cur]
        push!(path, cur)
    end
    reverse!(path)
    return (length(path) - 1, path)
end

"""
    connected_components(lat::AbstractLattice) -> Vector{Vector{Int}}

Return the connected components of the lattice graph, computed by BFS
over `neighbors(lat, ·)`. Each inner vector is sorted in ascending
order of site index, and components are returned in order of their
smallest member.

Useful for analysing diluted / defective lattices and for sanity
checks on user-defined `AbstractLattice` subtypes.

# Example
```julia
# A 3x3 open square stripped of one site is still connected:
length(connected_components(lat)) == 1
```
"""
function connected_components(lat::AbstractLattice)
    is_finite(lat) || throw(ArgumentError("connected_components requires a finite lattice"))
    N = num_sites(lat)
    visited = falses(N)
    comps = Vector{Vector{Int}}()
    for s in 1:N
        visited[s] && continue
        comp = Int[]
        queue = Int[s]
        visited[s] = true
        head = 1
        while head <= length(queue)
            u = queue[head]
            head += 1
            push!(comp, u)
            for v in neighbors(lat, u)
                if !visited[v]
                    visited[v] = true
                    push!(queue, v)
                end
            end
        end
        sort!(comp)
        push!(comps, comp)
    end
    return comps
end
