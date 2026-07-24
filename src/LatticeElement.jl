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

# ---- Generic element-center accessors --------------------------------
# Function declarations only; the actual methods need both
# AbstractLattice and Bond and live in `Bond.jl` (loaded after both).

"""
    num_elements(lat, e::AbstractLatticeElement) → Int

Number of geometric elements of centring `e` on `lat`. See
`Bond.jl` for the default `VertexCenter` / `BondCenter` methods.
"""
function num_elements end

"""
    elements(lat, e::AbstractLatticeElement)

Iterator over the underlying elements of centring `e` on `lat`. See
`Bond.jl` for the default `VertexCenter` / `BondCenter` methods.
"""
function elements end

"""
    element_position(lat, e::AbstractLatticeElement, i::Int)

Real-space position of the `i`-th element of centring `e` on `lat`.
See `Bond.jl` for the default `VertexCenter` / `BondCenter` methods.

# Indexing convention

The integer `i` follows the enumeration order of:

- `1:num_sites(lat)` for `VertexCenter`,
- `enumerate(bonds(lat))` for `BondCenter`,
- `enumerate(plaquettes(lat))` for `PlaquetteCenter`.

This means `i` is **lattice-specific**: a concrete lattice that
overrides `bonds(lat)` / `plaquettes(lat)` (or returns them in a
different order from another lattice type) implicitly redefines what
`i` means here. Generic code that round-trips via integer indices is
safe only within a single lattice instance; for cross-lattice
identification of an element, materialise it through `elements(lat, e)`
instead. A typed `BondIndex` / `PlaquetteIndex` wrapper that makes
this contract type-level is tracked as a follow-up.
"""
function element_position end

"""
    element_positions(lat, e::AbstractLatticeElement)

Iterator over the real-space positions of every element of centring
`e` on `lat`. See `Bond.jl` for the default implementation.
"""
function element_positions end
