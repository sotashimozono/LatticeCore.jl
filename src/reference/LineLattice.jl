"""
    LineLattice{T, B <: LatticeBoundary}(N, boundary)

A 1D linear chain of `N` sites with unit spacing. `boundary` is a
[`LatticeBoundary`](@ref) whose single axis component is one of
[`PeriodicAxis`](@ref), [`OpenAxis`](@ref), or [`TwistedAxis`](@ref).

LatticeCore's simplest reference implementation; paired with
[`SimpleSquareLattice`](@ref) it is also the canonical mock used by
downstream Monte Carlo unit tests.

# Examples
```julia
# Default: 1D periodic
line = LineLattice(5)

# Open boundary chain
chain = LineLattice(5, OpenAxis())

# Explicit LatticeBoundary (e.g. with a twist)
twisted = LineLattice(5, LatticeBoundary((TwistedAxis(π/4),)))
```
"""
struct LineLattice{T<:AbstractFloat,B<:LatticeBoundary} <: AbstractLattice{1,T}
    N::Int
    boundary::B
end

# Convenience constructors: default to Float64 positions. The axis BC
# is wrapped in a `LatticeBoundary` automatically so users can write
# `LineLattice(5, PeriodicAxis())` instead of the explicit tuple form.
function LineLattice(N::Int, axis::AbstractAxisBC=PeriodicAxis())
    return LineLattice(N, LatticeBoundary((axis,), NoModifier()))
end

function LineLattice(N::Int, boundary::B) where {B<:LatticeBoundary}
    return LineLattice{Float64,B}(N, boundary)
end

# ---- Required interface ----

num_sites(l::LineLattice) = l.N

position(l::LineLattice{T}, i::Int) where {T} = SVector{1,T}(T(i))

function neighbors(l::LineLattice, i::Int)
    ax = l.boundary.axes[1]
    ns = Int[]
    seen = Set{Int}()
    for step in (i + 1, i - 1)
        j, ok = apply_axis_bc(ax, step, l.N)
        if ok && j != i && !(j in seen)
            push!(ns, j)
            push!(seen, j)
        end
    end
    return ns
end

boundary(l::LineLattice) = l.boundary

size_trait(l::LineLattice) = FiniteSize((l.N,))

# ---- Bond iteration with wrapped (unit) displacement vectors ---------

function neighbor_bonds(l::LineLattice{T}, i::Int) where {T}
    ax = l.boundary.axes[1]
    out = Bond{1,T}[]
    seen = Set{Int}()
    for (step, dx) in ((i + 1, T(1)), (i - 1, T(-1)))
        j, ok = apply_axis_bc(ax, step, l.N)
        if ok && j != i && !(j in seen)
            push!(out, Bond{1,T}(i, j, SVector{1,T}(dx), :nearest))
            push!(seen, j)
        end
    end
    return out
end

function bonds(l::LineLattice{T}) where {T}
    return (b for i in 1:num_sites(l) for b in neighbor_bonds(l, i) if b.j > b.i)
end

# ---- Trait overrides ----

topology(::LineLattice) = TopologyTrait{:line}()

function periodicity(l::LineLattice)
    ax = l.boundary.axes[1]
    return ax isa OpenAxis ? Aperiodic() : Periodic()
end

# A 1D chain with PBC/TwistedBC is bipartite iff the cycle length is even.
# Open chain is always bipartite.
function is_bipartite(l::LineLattice)
    ax = l.boundary.axes[1]
    return ax isa OpenAxis ? true : iseven(l.N)
end

function reciprocal_support(l::LineLattice)
    ax = l.boundary.axes[1]
    return ax isa OpenAxis ? NoReciprocal() : HasReciprocal()
end

# ---- Coordinate conversions ----

"""
    to_real(lat::LineLattice, coord::LatticeCoord{1}) → RealSpace

Interpret the lattice coordinate as a unit-spacing Cartesian position.
"""
function to_real(::LineLattice{T}, coord::LatticeCoord{1}) where {T}
    return RealSpace{1,T}(SVector{1,T}(T(coord.cell[1])))
end

"""
    to_lattice(lat::LineLattice, rs::RealSpace{1}) → LatticeCoord{1}

Inverse of `to_real`: round the real-space x coordinate to the
nearest integer cell index.
"""
function to_lattice(::LineLattice, rs::RealSpace{1})
    return LatticeCoord((round(Int, rs.x[1]),), 1)
end
