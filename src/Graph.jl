"""
    Graph operations for `AbstractLattice`.

This file provides a small, dependency-light graph layer on top of the
`AbstractLattice` interface:

- [`adjacency_matrix`](@ref) — N × N adjacency as a `SparseMatrixCSC{Bool, Int}`
  (or a dense `Matrix{Bool}` when `sparse=false`),
- [`shortest_path`](@ref) — shortest path between two sites; unweighted
  BFS by default, weighted Dijkstra when a `weights` callback is given,
- [`distance_matrix`](@ref) — N × N all-pairs shortest-path distance
  matrix; Floyd–Warshall for small lattices, repeated Dijkstra otherwise,
- [`connected_components`](@ref) — connected components of the lattice
  graph as a `Vector{Vector{Int}}`,
- [`identity_weight`](@ref) — the default `weights` callback returning
  `1.0` for every bond.

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
    identity_weight(lat::AbstractLattice, bond::Bond) → 1.0

Default `weights` callback used by [`shortest_path`](@ref) and
[`distance_matrix`](@ref). Returns `1.0` for every bond.

This sentinel value is checked by `===` inside `shortest_path`:
when `weights === identity_weight`, the unweighted BFS code path is
selected for speed and to preserve the `Int`-distance return type
introduced in PR #44. Any other callback runs the Dijkstra code path
and yields `Float64` costs.

The callback signature is `(lat, bond) -> Real`. A custom weight
function should preserve this shape; for example, to use the
boundary-modifier [`bond_weight`](@ref) from
[`BoundaryCondition`](@ref AbstractBoundaryCondition):

```julia
ssd_w(lat, b) = bond_weight(boundary(lat).modifier, lat, b.i, b.j)
shortest_path(lat, src, dst; weights=ssd_w)
```

`identity_weight` itself ignores its arguments.
"""
identity_weight(::AbstractLattice, ::Bond) = 1.0

"""
    shortest_path(lat::AbstractLattice, src::Int, dst::Int;
                  weights::Function=identity_weight)
        -> (cost, path::Vector{Int})

Shortest path on the lattice graph from `src` to `dst`.

The implementation dispatches internally on the identity of the
`weights` callback:

- `weights === identity_weight` (the default) — the unweighted hop-
  count BFS introduced in PR #44 is used. Returns
  `(dist::Int, path::Vector{Int})`. `dist` is the number of bonds
  between `src` and `dst` (`0` if `src == dst`); on disconnect it
  returns `(typemax(Int), Int[])`. The no-keyword call form
  `shortest_path(lat, src, dst)` resolves to this branch and is fully
  backward compatible with PR #44.
- Any other callback — Dijkstra on the weighted graph. `weights`
  must be a callable `(lat, bond) -> Real` returning a non-negative
  edge cost. Returns `(cost::Float64, path::Vector{Int})`; on
  disconnect it returns `(Inf, Int[])`.

The `weights` callback is invoked once per directed edge traversal on
bonds produced by `neighbor_bonds(lat, u)`. The signature deliberately
takes a `Bond`, not just `(i, j)`, so callbacks can read `bond.type`,
`bond.vector`, etc. and so the same callback can be passed to
[`distance_matrix`](@ref). To re-use the existing
[`bond_weight`](@ref) of an [`AbstractBoundaryModifier`](@ref), wrap
it in a closure:

```julia
ssd_w(lat, b) = bond_weight(boundary(lat).modifier, lat, b.i, b.j)
shortest_path(lat, 1, 9; weights=ssd_w)
```

Negative-weight edges are not supported; behaviour is undefined if
`weights` returns a negative value.

# Examples
```julia
lat = SimpleSquareLattice(3, 3, OpenAxis())
d, p = shortest_path(lat, 1, 9)                       # BFS, Int distance
c, q = shortest_path(lat, 1, 9; weights=(lat, b) -> 1.0)  # Dijkstra, Float64 cost
c == Float64(d)                                        # true
```
"""
function shortest_path(
    lat::AbstractLattice, src::Int, dst::Int; weights::Function=identity_weight
)
    # Method dispatch ignores keyword arguments, so we can only have
    # one `shortest_path(lat, src, dst)` method. We branch internally on
    # the `weights` callback:
    #
    # - `weights === identity_weight` (default) keeps the PR #44 BFS
    #   call shape: returns `(dist::Int, path::Vector{Int})` and
    #   `(typemax(Int), Int[])` on disconnect.
    # - Any other callback runs Dijkstra and returns
    #   `(cost::Float64, path::Vector{Int})` with `(Inf, Int[])` on
    #   disconnect.
    #
    # This preserves backward compatibility with the unweighted call
    # while letting users opt into Float64-cost Dijkstra by passing a
    # custom `weights`.
    if weights === identity_weight
        return _shortest_path_bfs(lat, src, dst)
    else
        return _shortest_path_dijkstra(lat, src, dst, weights)
    end
end

# Internal: BFS implementation. Kept separate from the Dijkstra path so
# the no-`weights` call shape (PR #44 API) keeps its `Int` distance
# return type and avoids paying the priority-queue overhead.
function _shortest_path_bfs(lat::AbstractLattice, src::Int, dst::Int)
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

# Internal: Dijkstra single-source, single-target. Uses a binary heap
# (O((V+E) log V)) on `(cost, site)` pairs. We index the heap by site
# count which is small for any lattice we expect; if we ever need
# decrease-key, we can switch to an indexed heap. The "stale entry"
# trick (skip entries whose recorded cost no longer matches `dist`)
# is the standard workaround for a non-indexed binary heap.
function _shortest_path_dijkstra(
    lat::AbstractLattice, src::Int, dst::Int, weights::Function
)
    is_finite(lat) || throw(ArgumentError("shortest_path requires a finite lattice"))
    N = num_sites(lat)
    (1 <= src <= N) || throw(BoundsError(lat, src))
    (1 <= dst <= N) || throw(BoundsError(lat, dst))

    if src == dst
        return (0.0, [src])
    end

    dist, parent = _dijkstra(lat, src, weights, dst)

    if dist[dst] == Inf
        return (Inf, Int[])
    end

    path = Int[dst]
    cur = dst
    while cur != src
        cur = parent[cur]
        push!(path, cur)
    end
    reverse!(path)
    return (dist[dst], path)
end

# Internal: single-source Dijkstra over `lat`, using `weights(lat, b)`
# for each `Bond` returned by `neighbor_bonds(lat, u)`. Returns the
# full `(dist, parent)` arrays so callers can either ask for a single
# `dst` (early termination via the optional `target` argument) or build
# an all-pairs distance matrix on top.
function _dijkstra(lat::AbstractLattice, src::Int, weights::Function, target::Int=0)
    N = num_sites(lat)
    dist = fill(Inf, N)
    parent = fill(0, N)
    dist[src] = 0.0

    # Heap stores `(cost, site)` tuples. Default `Base.isless` on
    # tuples is lex-order, so `cost` dominates which is what we want.
    heap = Tuple{Float64,Int}[(0.0, src)]
    while !isempty(heap)
        d, u = _heappop!(heap)
        # Stale entry guard: another, cheaper path to `u` has been
        # popped already, so this entry is obsolete.
        d > dist[u] && continue
        u == target && break
        for b in neighbor_bonds(lat, u)
            v = b.j
            w = Float64(weights(lat, b))
            alt = d + w
            if alt < dist[v]
                dist[v] = alt
                parent[v] = u
                _heappush!(heap, (alt, v))
            end
        end
    end
    return dist, parent
end

# Internal: minimal binary-heap push/pop on a Vector. We keep this
# self-contained (no DataStructures.jl dep) because Graph.jl is meant
# to be dependency-light. `_heappush!` / `_heappop!` operate on the
# default `Base.isless` ordering.
function _heappush!(heap::Vector{T}, x::T) where {T}
    push!(heap, x)
    i = length(heap)
    while i > 1
        parent = i >> 1
        if isless(heap[i], heap[parent])
            heap[i], heap[parent] = heap[parent], heap[i]
            i = parent
        else
            break
        end
    end
    return heap
end

function _heappop!(heap::Vector{T}) where {T}
    top = heap[1]
    last = pop!(heap)
    if !isempty(heap)
        heap[1] = last
        n = length(heap)
        i = 1
        while true
            l = 2i
            r = 2i + 1
            smallest = i
            if l <= n && isless(heap[l], heap[smallest])
                smallest = l
            end
            if r <= n && isless(heap[r], heap[smallest])
                smallest = r
            end
            smallest == i && break
            heap[i], heap[smallest] = heap[smallest], heap[i]
            i = smallest
        end
    end
    return top
end

"""
    distance_matrix(lat::AbstractLattice;
                    weights::Function=identity_weight) → Matrix{Float64}

All-pairs shortest-path distance matrix on the lattice graph. Returns
an `N × N` `Matrix{Float64}` with `D[i, j]` equal to the cost of the
shortest path from site `i` to site `j` under `weights`. Diagonal
entries are `0.0`. Pairs of sites in distinct connected components
are filled with `Inf`.

The default `weights=identity_weight` reproduces unweighted hop-count
distances (`Float64`-promoted for uniformity with the weighted case).
The callback signature is `(lat, bond) -> Real` — see
[`identity_weight`](@ref) for the convention and an SSD example.

# Algorithm

For lattices with `num_sites(lat) <= floyd_warshall_threshold` (default
`64`) the implementation uses dense Floyd–Warshall (`O(N^3)`,
arithmetic-only). For larger lattices it runs Dijkstra from each
source (`O(N · (V + E) log V)`), which is asymptotically cheaper on
sparse lattice graphs.

# Keyword arguments

- `weights::Function=identity_weight` — edge cost callback as above.
- `floyd_warshall_threshold::Int=64` — switch-over `N` between the
  two algorithms. Mostly a tuning knob; the default is conservative.

# Example

```julia
lat = SimpleSquareLattice(4, 4, OpenAxis())
D = distance_matrix(lat)
@assert D == transpose(D)            # symmetric for undirected graph
@assert all(iszero, diag(D))
```
"""
function distance_matrix(
    lat::AbstractLattice;
    weights::Function=identity_weight,
    floyd_warshall_threshold::Int=64,
)
    is_finite(lat) || throw(ArgumentError("distance_matrix requires a finite lattice"))
    N = num_sites(lat)
    if N <= floyd_warshall_threshold
        return _distance_matrix_floyd_warshall(lat, weights, N)
    else
        return _distance_matrix_dijkstra(lat, weights, N)
    end
end

# Internal: dense Floyd–Warshall. We seed the matrix from the bond
# iterator directly (rather than going through `adjacency_matrix`)
# because we need real-valued weights, not the Bool adjacency.
# Multi-edges are folded by taking the minimum (consistent with
# `adjacency_matrix`'s `combine = |` for Bool).
function _distance_matrix_floyd_warshall(lat::AbstractLattice, weights::Function, N::Int)
    D = fill(Inf, N, N)
    @inbounds for i in 1:N
        D[i, i] = 0.0
    end
    for b in elements(lat, BondCenter())
        i, j = b.i, b.j
        i == j && continue
        w = Float64(weights(lat, b))
        # Symmetric: bonds are undirected. If a lattice publishes both
        # orientations (i->j and j->i) with different weights, the
        # smaller one wins on each side. We do not currently encode
        # directed bonds.
        if w < D[i, j]
            D[i, j] = w
            D[j, i] = w
        end
    end
    @inbounds for k in 1:N, i in 1:N
        dik = D[i, k]
        dik == Inf && continue
        for j in 1:N
            dkj = D[k, j]
            dkj == Inf && continue
            alt = dik + dkj
            if alt < D[i, j]
                D[i, j] = alt
            end
        end
    end
    return D
end

# Internal: repeated single-source Dijkstra. Faster than Floyd–Warshall
# on sparse lattice graphs once N is large (graphs from
# `AbstractLattice` always have O(N) edges).
function _distance_matrix_dijkstra(lat::AbstractLattice, weights::Function, N::Int)
    D = Matrix{Float64}(undef, N, N)
    for s in 1:N
        dist, _ = _dijkstra(lat, s, weights)
        @inbounds for j in 1:N
            D[s, j] = dist[j]
        end
    end
    return D
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
