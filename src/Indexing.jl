"""
    AbstractIndexing

Abstract supertype for strategies that linearise a
[`LatticeCoord`](@ref) into a 1-based site index. Concrete subtypes
determine the ordering of sites along axes and the placement of
sublattices within each cell.

The indexing method is deliberately **decoupled** from the coordinate
system: a lattice may switch its `AbstractIndexing` without changing
how it describes positions. See
`dev/note/04_architecture/03_boundary_and_coordinates/README.md` for
the design rationale.

# Required interface (per concrete indexing / dimension)

- `site_index(indexing, dims::NTuple{D, Int}, nsub::Int, coord::LatticeCoord{D})::Int`
- `lattice_coord(indexing, dims::NTuple{D, Int}, nsub::Int, i::Int)::LatticeCoord{D}`

These two methods must round-trip: `site_index ∘ lattice_coord` is
the identity on `1:prod(dims) * nsub`.

# Sublattice convention

The sublattice is the **innermost** (fastest) index: sites within a
single cell are contiguous. For `nsub == 1` this reduces to the usual
single-sublattice formulas and can be read at a glance.
"""
abstract type AbstractIndexing end

"""Row-major (C-style): x is the fastest axis within a row."""
struct RowMajor <: AbstractIndexing end

"""Column-major (Fortran-style): y is the fastest axis."""
struct ColMajor <: AbstractIndexing end

"""
    Snake

Boustrophedon (snake) indexing: row-major layout with alternating
row direction. Useful when physical adjacency between neighbouring
indices matters (e.g. for certain serialisation layouts).
"""
struct Snake <: AbstractIndexing end

# ---- Generic functions ----------------------------------------------

"""
    site_index(indexing, dims, nsub, coord) → Int

Linearise `coord` into a 1-based site index. See [`AbstractIndexing`](@ref)
for the required contract.
"""
function site_index end

"""
    lattice_coord(indexing, dims, nsub, i) → LatticeCoord

Inverse of [`site_index`](@ref).
"""
function lattice_coord end

# ---- 1D ---------------------------------------------------------------
#
# In 1D there is no choice of ordering: RowMajor / ColMajor / Snake
# all reduce to "walk along the single axis in order", so we define
# only RowMajor for 1D and let the others fall through (they would
# simply error out; 1D Snake is actively meaningless).

function site_index(::RowMajor, dims::NTuple{1,Int}, nsub::Int, coord::LatticeCoord{1})
    cx = coord.cell[1]
    s = coord.sublattice
    return (cx - 1) * nsub + s
end

function lattice_coord(::RowMajor, ::NTuple{1,Int}, nsub::Int, i::Int)
    c = (i - 1) ÷ nsub
    s = (i - 1) % nsub + 1
    return LatticeCoord((c + 1,), s)
end

# ---- 2D RowMajor ------------------------------------------------------
#
# site_index(cx, cy, s) = ((cy - 1) * Lx + (cx - 1)) * nsub + s

function site_index(::RowMajor, dims::NTuple{2,Int}, nsub::Int, coord::LatticeCoord{2})
    cx, cy = coord.cell
    s = coord.sublattice
    Lx = dims[1]
    return ((cy - 1) * Lx + (cx - 1)) * nsub + s
end

function lattice_coord(::RowMajor, dims::NTuple{2,Int}, nsub::Int, i::Int)
    c = (i - 1) ÷ nsub
    s = (i - 1) % nsub + 1
    Lx = dims[1]
    cx = c % Lx + 1
    cy = c ÷ Lx + 1
    return LatticeCoord((cx, cy), s)
end

# ---- 2D ColMajor ------------------------------------------------------
#
# site_index(cx, cy, s) = ((cx - 1) * Ly + (cy - 1)) * nsub + s

function site_index(::ColMajor, dims::NTuple{2,Int}, nsub::Int, coord::LatticeCoord{2})
    cx, cy = coord.cell
    s = coord.sublattice
    Ly = dims[2]
    return ((cx - 1) * Ly + (cy - 1)) * nsub + s
end

function lattice_coord(::ColMajor, dims::NTuple{2,Int}, nsub::Int, i::Int)
    c = (i - 1) ÷ nsub
    s = (i - 1) % nsub + 1
    Ly = dims[2]
    cy = c % Ly + 1
    cx = c ÷ Ly + 1
    return LatticeCoord((cx, cy), s)
end

# ---- 2D Snake ---------------------------------------------------------
#
# Row-major layout, but within each row the x-direction alternates:
# odd rows walk x = 1..Lx, even rows walk x = Lx..1.

function site_index(::Snake, dims::NTuple{2,Int}, nsub::Int, coord::LatticeCoord{2})
    cx, cy = coord.cell
    s = coord.sublattice
    Lx = dims[1]
    effective_x = isodd(cy) ? cx : (Lx + 1 - cx)
    cell_index = (cy - 1) * Lx + effective_x
    return (cell_index - 1) * nsub + s
end

function lattice_coord(::Snake, dims::NTuple{2,Int}, nsub::Int, i::Int)
    c = (i - 1) ÷ nsub
    s = (i - 1) % nsub + 1
    Lx = dims[1]
    cell_index = c + 1
    cy = (cell_index - 1) ÷ Lx + 1
    effective_x = (cell_index - 1) % Lx + 1
    cx = isodd(cy) ? effective_x : (Lx + 1 - effective_x)
    return LatticeCoord((cx, cy), s)
end
