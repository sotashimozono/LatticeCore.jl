"""
    SimpleSquareLattice{T, B <: LatticeBoundary}(Lx, Ly, boundary)

A 2D square lattice of `Lx × Ly` sites with unit spacing. `boundary`
is a [`LatticeBoundary`](@ref) whose two axis components independently
select between [`PeriodicAxis`](@ref), [`OpenAxis`](@ref), and
[`TwistedAxis`](@ref) — so mixed BCs (e.g. cylinders) are supported
natively.

LatticeCore's 2D reference implementation. It exists so that the core
interface can be exercised end-to-end without depending on
`Lattice2D.jl`. The production-grade 2D lattice, with the full
topology catalogue (triangular, honeycomb, kagome, ...), sublattice
site types, and reciprocal-lattice machinery, lives in `Lattice2D.jl`.

# Site indexing

Sites are laid out in row-major order:

    site_index(x, y) = (y - 1) * Lx + x     # 1 <= x <= Lx, 1 <= y <= Ly

so sites 1..Lx form the bottom row, Lx+1..2Lx the next, and so on.

# Examples
```julia
# Default: 2D PBC
sq = SimpleSquareLattice(3, 3)

# Uniform open boundary
open_sq = SimpleSquareLattice(3, 3, OpenAxis())

# Cylinder: periodic in x, open in y
cylinder = SimpleSquareLattice(3, 3, LatticeBoundary((PeriodicAxis(), OpenAxis())))
```
"""
struct SimpleSquareLattice{T<:AbstractFloat,B<:LatticeBoundary,L<:AbstractSiteLayout} <:
       AbstractLattice{2,T}
    Lx::Int
    Ly::Int
    boundary::B
    layout::L
end

# Convenience constructors. The single-axis overload applies the same
# BC to both axes; mixed cases go through `LatticeBoundary` directly.
# The site layout defaults to a uniform Ising layout.
function SimpleSquareLattice(
    Lx::Int,
    Ly::Int,
    axis::AbstractAxisBC=PeriodicAxis();
    layout::AbstractSiteLayout=UniformLayout(IsingSite()),
)
    return SimpleSquareLattice(Lx, Ly, LatticeBoundary((axis, axis), NoModifier()); layout)
end

function SimpleSquareLattice(
    Lx::Int,
    Ly::Int,
    boundary::LatticeBoundary;
    layout::AbstractSiteLayout=UniformLayout(IsingSite()),
)
    return SimpleSquareLattice{Float64,typeof(boundary),typeof(layout)}(
        Lx, Ly, boundary, layout
    )
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

# Single helper shared with `neighbor_bonds`: walk the four axis steps,
# apply per-axis BCs, deduplicate.
function _square_steps(l::SimpleSquareLattice, i::Int)
    x, y = _site_to_xy(l, i)
    bx, by = l.boundary.axes
    steps = Tuple{Int,Int,Int,Int}[]           # (dx, dy, nx, ny)
    for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
        nx_raw = x + dx
        ny_raw = y + dy
        nx, ok_x = apply_axis_bc(bx, nx_raw, l.Lx)
        ok_x || continue
        ny, ok_y = apply_axis_bc(by, ny_raw, l.Ly)
        ok_y || continue
        push!(steps, (dx, dy, nx, ny))
    end
    return steps
end

function neighbors(l::SimpleSquareLattice, i::Int)
    ns = Int[]
    seen = Set{Int}()
    for (_, _, nx, ny) in _square_steps(l, i)
        j = _xy_to_site(l, nx, ny)
        if j != i && !(j in seen)
            push!(ns, j)
            push!(seen, j)
        end
    end
    return ns
end

boundary(l::SimpleSquareLattice) = l.boundary

site_layout(l::SimpleSquareLattice) = l.layout

size_trait(l::SimpleSquareLattice) = FiniteSize((l.Lx, l.Ly))

# ---- Bond iteration with wrapped (unit) displacement vectors ---------

function neighbor_bonds(l::SimpleSquareLattice{T}, i::Int) where {T}
    out = Bond{2,T}[]
    seen = Set{Int}()
    for (dx, dy, nx, ny) in _square_steps(l, i)
        j = _xy_to_site(l, nx, ny)
        if j != i && !(j in seen)
            vec = SVector{2,T}(T(dx), T(dy))
            push!(out, Bond{2,T}(i, j, vec, :nearest))
            push!(seen, j)
        end
    end
    return out
end

function bonds(l::SimpleSquareLattice{T}) where {T}
    return (b for i in 1:num_sites(l) for b in neighbor_bonds(l, i) if b.j > b.i)
end

# ---- Trait overrides ----

topology(::SimpleSquareLattice) = TopologyTrait{:square}()

function periodicity(l::SimpleSquareLattice)
    all_periodic = all(!(a isa OpenAxis) for a in l.boundary.axes)
    return all_periodic ? Periodic() : Aperiodic()
end

# Square lattice graph bipartiteness:
# - OBC axes never introduce odd cycles.
# - PBC / Twisted cycles must have even length for bipartiteness.
# We require every axis to be bipartite in isolation; the combined
# graph is bipartite iff every axis along which the lattice wraps has
# even length.
function is_bipartite(l::SimpleSquareLattice)
    bx, by = l.boundary.axes
    bx_ok = (bx isa OpenAxis) || iseven(l.Lx)
    by_ok = (by isa OpenAxis) || iseven(l.Ly)
    return bx_ok && by_ok
end

function reciprocal_support(l::SimpleSquareLattice)
    all_periodic = all(!(a isa OpenAxis) for a in l.boundary.axes)
    return all_periodic ? HasReciprocal() : NoReciprocal()
end

# ---- Coordinate conversions ----

"""
    to_real(lat::SimpleSquareLattice, coord::LatticeCoord{2}) → RealSpace

Interpret the lattice coordinate as a unit-spacing Cartesian position.
"""
function to_real(::SimpleSquareLattice{T}, coord::LatticeCoord{2}) where {T}
    return RealSpace{2,T}(SVector{2,T}(T(coord.cell[1]), T(coord.cell[2])))
end

"""
    to_lattice(lat::SimpleSquareLattice, rs::RealSpace{2}) → LatticeCoord{2}

Inverse of `to_real`: round each real-space component to the nearest
integer cell index.
"""
function to_lattice(::SimpleSquareLattice, rs::RealSpace{2})
    return LatticeCoord((round(Int, rs.x[1]), round(Int, rs.x[2])), 1)
end

# ---- Reciprocal lattice ----

"""
    basis_vectors(lat::SimpleSquareLattice) → SMatrix{2, 2, T}

Real-space basis matrix. For the unit-spacing reference square this
is the 2×2 identity.
"""
function basis_vectors(::SimpleSquareLattice{T}) where {T}
    return SMatrix{2,2,T}(one(T), zero(T), zero(T), one(T))
end

# ---- Translation-cell motif ------------------------------------------
#
# The finite square lattice shares its unit-cell motif with
# `InfiniteSquareLattice`, so a TN builder written against
# `site_orbits` / `bond_orbits` runs unchanged on both. Finiteness and
# the boundary condition only affect how the motif is tiled and wrapped,
# not the motif itself.

translation_vectors(l::SimpleSquareLattice) = basis_vectors(l)

num_basis_sites(::SimpleSquareLattice) = 1

function cell_bonds(::SimpleSquareLattice)
    return (CellBond(1, 1, (1, 0), :nearest), CellBond(1, 1, (0, 1), :nearest))
end

function reciprocal_lattice(lat::SimpleSquareLattice{T}) where {T}
    reciprocal_support(lat) isa HasReciprocal || throw(
        ArgumentError(
            "SimpleSquareLattice has no reciprocal lattice unless both axes are periodic",
        ),
    )
    A = basis_vectors(lat)
    B = SMatrix{2,2,T}(T(2π) * inv(transpose(A)))
    return monkhorst_pack(B, (lat.Lx, lat.Ly))
end

# ---- FFT fast-path opt-in (see `src/StructureFactor.jl`) -------------
#
# `SimpleSquareLattice` uses row-major site indexing:
# `site_index(x, y) = (y - 1) * Lx + x`. Julia's column-major reshape
# over `(Lx, Ly)` yields `M[x, y] = state[(y - 1) * Lx + x]`, which
# lines up with the natural grid (x along axis 1, y along axis 2).
# Registered here so `Lattice2D` / downstream code can rely on the
# opt-in without depending on FFTW.
_has_known_grid(::SimpleSquareLattice) = true
_reshape_state(::SimpleSquareLattice, state, dims) = reshape(state, dims)
