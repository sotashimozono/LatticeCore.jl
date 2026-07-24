"""
    InfiniteSquareLattice{T, L <: AbstractSiteLayout}(; layout)

A truly infinite 2D square lattice with unit spacing — the
thermodynamic-limit counterpart of [`SimpleSquareLattice`](@ref) under
periodic boundary conditions.

It carries no size: `size_trait` is [`InfiniteSize`](@ref) and
[`num_sites`](@ref) throws. It is *not* meant to be walked with the
linear `1:num_sites` site API. Instead it is described by its unit-cell
motif and accessed lazily:

- [`site_orbits`](@ref) / [`bond_orbits`](@ref) — the finite
  fundamental domain (one basis site, two bonds `+x` and `+y`);
- [`cell_position`](@ref), [`neighbors_at`](@ref),
  [`incident_cell_bonds`](@ref) — on-demand access to any concrete
  [`CellSite`](@ref), computed without materialising the lattice.

This is the substrate for building an infinite tensor network (e.g. the
Ising partition function in the thermodynamic limit): place one tensor
per site orbit and one per bond orbit, then contract along the motif.

# Bridge to a finite lattice

If a finite approximation is ever needed, [`materialize`](@ref) tiles
the motif into a periodic [`SimpleSquareLattice`](@ref):

```julia
inf = InfiniteSquareLattice()
fin = materialize(inf; dims = (8, 8))   # 8×8 PBC SimpleSquareLattice
```

but the lazy accessors above never require this step.

# Examples
```julia
inf = InfiniteSquareLattice()
site_orbits(inf)                       # 1:1  (one basis site)
collect(bond_orbits(inf))              # +x and +y CellBonds
s = CellSite((3, -2))                  # cell (3, -2), basis 1
cell_position(inf, s)                  # SVector(3.0, -2.0)
neighbors_at(inf, s)                   # the four nearest neighbours
```
"""
struct InfiniteSquareLattice{T<:AbstractFloat,L<:AbstractSiteLayout} <: AbstractLattice{2,T}
    layout::L
end

function InfiniteSquareLattice(; layout::AbstractSiteLayout=UniformLayout(IsingSite()))
    return InfiniteSquareLattice{Float64,typeof(layout)}(layout)
end

# ---- Translation-cell interface --------------------------------------

function translation_vectors(::InfiniteSquareLattice{T}) where {T}
    return SMatrix{2,2,T}(one(T), zero(T), zero(T), one(T))
end

num_basis_sites(::InfiniteSquareLattice) = 1

function cell_bonds(::InfiniteSquareLattice)
    return (CellBond(1, 1, (1, 0), :nearest), CellBond(1, 1, (0, 1), :nearest))
end

# ---- Size / traits ---------------------------------------------------

size_trait(::InfiniteSquareLattice) = InfiniteSize()

"""
    num_sites(::InfiniteSquareLattice)

Always throws `DomainError`: an infinite lattice has no finite site
count. Use [`site_orbits`](@ref) for the fundamental domain, or
[`materialize`](@ref) to obtain a finite sample.
"""
function num_sites(::InfiniteSquareLattice)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteSquareLattice has no finite site count; use `site_orbits` " *
            "for the unit-cell basis or `materialize(lat; dims)` for a finite sample",
        ),
    )
end

topology(::InfiniteSquareLattice) = TopologyTrait{:square}()
periodicity(::InfiniteSquareLattice) = Periodic()
is_bipartite(::InfiniteSquareLattice) = true
reciprocal_support(::InfiniteSquareLattice) = HasReciprocal()
site_layout(l::InfiniteSquareLattice) = l.layout

# Conceptually the thermodynamic limit of full PBC.
boundary(::InfiniteSquareLattice) = LatticeBoundary((PeriodicAxis(), PeriodicAxis()))

"""
    basis_vectors(::InfiniteSquareLattice) → SMatrix{2, 2, T}

Real-space basis matrix (the 2×2 identity for unit spacing); equal to
[`translation_vectors`](@ref) for this single-basis lattice.
"""
basis_vectors(l::InfiniteSquareLattice) = translation_vectors(l)

# ---- Linear-index API is undefined for an infinite lattice -----------

function position(::InfiniteSquareLattice, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteSquareLattice has no linear site index; address sites by " *
            "`CellSite` and use `cell_position(lat, site)`",
        ),
    )
end

function neighbors(::InfiniteSquareLattice, ::Int)
    return throw(
        DomainError(
            InfiniteSize(),
            "InfiniteSquareLattice has no linear site index; use " *
            "`neighbors_at(lat, ::CellSite)`",
        ),
    )
end

# ---- Bridge: materialise into a finite periodic sample ---------------

"""
    materialize(lat::InfiniteSquareLattice; dims::NTuple{2, Int}) → SimpleSquareLattice

Tile the motif into a `dims[1] × dims[2]` periodic
[`SimpleSquareLattice`](@ref), carrying the same site layout. This is
the optional infinite-abstract → finite-materialisation bridge; the
lazy translation-cell accessors do not require it.
"""
function materialize(lat::InfiniteSquareLattice; dims::NTuple{2,Int})
    return SimpleSquareLattice(dims[1], dims[2], PeriodicAxis(); layout=lat.layout)
end
