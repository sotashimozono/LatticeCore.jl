"""
    PlaquetteRule

Declarative description of a plaquette shape anchored in a unit cell,
used as the input to the per-sample plaquette enumeration. A concrete
topology stores one `PlaquetteRule` per plaquette kind (e.g. Square has
one rule, Triangular has two — up-triangle and down-triangle, Kagome
has three — up-triangle, down-triangle, hexagon).

# Fields
- `corners::Vector{NTuple{3, Int}}` — the boundary vertices of the
  plaquette in **cyclic order**, each expressed as
  `(sublattice_id, dx, dy)`. `sublattice_id` is the 1-based sublattice
  index within the unit cell; `(dx, dy)` is the cell offset relative
  to the anchor cell. Every rule must include at least one corner at
  `(_, 0, 0)` (the anchor).
- `type::Symbol` — bond-type-like tag for downstream dispatch
  (`:square`, `:up_triangle`, `:down_triangle`, `:hexagon`, …).

# Example

```julia
# Unit square on a single-sublattice square lattice:
square_rule = PlaquetteRule(
    [(1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)],
    :square,
)
```

Mirrors the `Connection` / `Bond` pattern: `PlaquetteRule` is the
topology-level template, `Plaquette` is the per-sample materialised
value produced by applying the rule under a boundary condition.
"""
struct PlaquetteRule
    corners::Vector{NTuple{3,Int}}
    type::Symbol
end

"""
    Plaquette{D, T}

A plaquette (face) on a concrete lattice, materialised from a
[`PlaquetteRule`](@ref).

# Fields
- `vertices::Vector{Int}` — 1-based site indices of the boundary
  vertices, in cyclic order
- `center::SVector{D, T}` — real-space centroid
- `type::Symbol` — tag inherited from the rule

Obtained via `plaquettes(lat)` / `neighbor_plaquettes(lat, i)` /
`element_position(lat, PlaquetteCenter(), p)` or the generic
element-center API.
"""
struct Plaquette{D,T}
    vertices::Vector{Int}
    center::SVector{D,T}
    type::Symbol
end

"""
    plaquettes(lat::AbstractLattice)

Iterator over all plaquettes on `lat`. The default implementation
throws `MethodError` — concrete lattices that have a notion of
plaquettes must implement this method (e.g. by walking cells ×
`PlaquetteRule`s under their boundary condition).
"""
function plaquettes end

"""
    neighbor_plaquettes(lat::AbstractLattice, i::Int)

Iterator over plaquettes that have site `i` on their boundary.
Default implementation filters `plaquettes(lat)` by membership; concrete
lattices may override for efficiency.
"""
function neighbor_plaquettes(lat::AbstractLattice, i::Int)
    return (p for p in plaquettes(lat) if i in p.vertices)
end

"""
    plaquette_center(p::Plaquette) → SVector

Geometric center of the plaquette. For a materialised `Plaquette`
this is just a field read.
"""
plaquette_center(p::Plaquette) = p.center

# ---- Element-center defaults for PlaquetteCenter --------------------

# Default: count via the plaquettes iterator. Concrete lattices may
# override with an O(1) cell×kind formula.
function num_elements(lat::AbstractLattice, ::PlaquetteCenter)
    return count(_ -> true, plaquettes(lat))
end

# Default: return the iterator as-is.
elements(lat::AbstractLattice, ::PlaquetteCenter) = plaquettes(lat)

# Default: materialise and index. O(num_plaquettes) per call; concrete
# lattices may override for O(1).
function element_position(lat::AbstractLattice, ::PlaquetteCenter, i::Int)
    ps = collect(plaquettes(lat))
    return ps[i].center
end

# ---- Same-centring adjacency (line graph / dual graph) --------------

"""
    element_neighbors(lat, e::AbstractLatticeElement, i::Int)

Neighbours of the `i`-th element of centring `e`, under the adjacency
appropriate to `e`:

- `VertexCenter` → `neighbors(lat, i)` (the usual graph)
- `BondCenter` → **line graph**: other bonds sharing a vertex with bond `i`
- `PlaquetteCenter` → **dual graph**: plaquettes sharing an edge with plaquette `i`
- `CellCenter` → throws unless overridden

Concrete lattices may override any of these for efficiency or to
support custom adjacency rules (e.g. "only same-type bonds").
"""
function element_neighbors end

element_neighbors(lat::AbstractLattice, ::VertexCenter, i::Int) = neighbors(lat, i)

function element_neighbors(lat::AbstractLattice, ::BondCenter, i::Int)
    bs = collect(bonds(lat))
    b = bs[i]
    endpoints = (b.i, b.j)
    out = Int[]
    for (k, other) in enumerate(bs)
        k == i && continue
        if other.i in endpoints || other.j in endpoints
            push!(out, k)
        end
    end
    return out
end

function element_neighbors(lat::AbstractLattice, ::PlaquetteCenter, i::Int)
    ps = collect(plaquettes(lat))
    p = ps[i]
    p_edges = _boundary_edges(p)
    out = Int[]
    for (k, other) in enumerate(ps)
        k == i && continue
        any_shared = false
        for e in p_edges
            if e in _boundary_edges(other) || reverse(e) in _boundary_edges(other)
                any_shared = true
                break
            end
        end
        any_shared && push!(out, k)
    end
    return out
end

# Boundary edges of a plaquette, as unordered `Tuple{Int,Int}` pairs
# walked cyclically around the vertex list.
function _boundary_edges(p::Plaquette)
    vs = p.vertices
    n = length(vs)
    return [(vs[i], vs[mod1(i + 1, n)]) for i in 1:n]
end

# ---- Cross-centring incidence ---------------------------------------

"""
    incident(lat, from::AbstractLatticeElement, to::AbstractLatticeElement, i::Int)

Return the integer indices of elements of centring `to` that are
incident to the `i`-th element of centring `from`. Defaults cover:

- `VertexCenter` ↔ `BondCenter`: bond-site incidence
- `VertexCenter` ↔ `PlaquetteCenter`: site-plaquette incidence
- `BondCenter` ↔ `PlaquetteCenter`: bond-plaquette incidence

Same-centring pairs (`from == to`) fall through to
[`element_neighbors`](@ref) so that `incident(lat, E(), E(), i)` is the
adjacency under centring `E`.

Concrete lattices may override any pair for O(1) access — the default
implementations here are O(num_elements) materialisations meant for
correctness, not hot-path use.
"""
function incident end

# Same-centring → adjacency shortcut.
function incident(
    lat::AbstractLattice, e::AbstractLatticeElement, e2::AbstractLatticeElement, i::Int
)
    e === e2 || error("no default for incident($(typeof(e)) → $(typeof(e2)))")
    return element_neighbors(lat, e, i)
end

# Vertex → Bond: bonds touching site i
function incident(lat::AbstractLattice, ::VertexCenter, ::BondCenter, i::Int)
    bs = collect(bonds(lat))
    return [k for (k, b) in enumerate(bs) if b.i == i || b.j == i]
end

# Bond → Vertex: endpoints of bond i
function incident(lat::AbstractLattice, ::BondCenter, ::VertexCenter, i::Int)
    b = collect(bonds(lat))[i]
    return [b.i, b.j]
end

# Vertex → Plaquette: plaquettes containing site i
function incident(lat::AbstractLattice, ::VertexCenter, ::PlaquetteCenter, i::Int)
    ps = collect(plaquettes(lat))
    return [k for (k, p) in enumerate(ps) if i in p.vertices]
end

# Plaquette → Vertex: boundary vertices of plaquette i (in cyclic order)
function incident(lat::AbstractLattice, ::PlaquetteCenter, ::VertexCenter, i::Int)
    return collect(plaquettes(lat))[i].vertices
end

# Bond → Plaquette: plaquettes whose boundary contains bond i
function incident(lat::AbstractLattice, ::BondCenter, ::PlaquetteCenter, i::Int)
    bs = collect(bonds(lat))
    b = bs[i]
    target = (b.i, b.j)
    ps = collect(plaquettes(lat))
    out = Int[]
    for (k, p) in enumerate(ps)
        for e in _boundary_edges(p)
            if e == target || e == reverse(target)
                push!(out, k)
                break
            end
        end
    end
    return out
end

# Plaquette → Bond: bonds forming the boundary of plaquette i, as
# integer bond indices. Each boundary edge is matched against
# `bonds(lat)` by endpoint set; unmatched edges are skipped.
function incident(lat::AbstractLattice, ::PlaquetteCenter, ::BondCenter, i::Int)
    p = collect(plaquettes(lat))[i]
    bs = collect(bonds(lat))
    out = Int[]
    for e in _boundary_edges(p)
        for (k, b) in enumerate(bs)
            if (b.i, b.j) == e || (b.j, b.i) == e
                push!(out, k)
                break
            end
        end
    end
    return out
end
