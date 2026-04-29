"""
    AbstractMomentumLattice{D, T} <: AbstractLattice{D, T}

Abstract supertype for k-space lattices. A momentum lattice is itself
an `AbstractLattice` — its "sites" are k-points, its "positions" are
k-vectors in the reciprocal basis — so 02's lattice interface
(`num_sites`, `position`, traits, test suite) can be re-used for
k-space work. See `dev/note/04_architecture/05_momentum_space`.

# Required interface for concrete subtypes
- `num_k_points(ml)::Int`
- `k_point(ml, i::Int)::SVector{D, T}`
- `reciprocal_basis(ml)::SMatrix{D, D, T}`

The graph-neighbour side of `AbstractLattice` (`neighbors`, `bonds`)
is deliberately vacuous for momentum lattices: k-space is not a graph
in the same sense. Concrete momentum lattices return an empty
neighbour list; MC code should never walk them as a graph.
"""
abstract type AbstractMomentumLattice{D,T} <: AbstractLattice{D,T} end

"""
    num_k_points(ml::AbstractMomentumLattice) → Int

Number of k-points stored in the lattice.
"""
function num_k_points end

"""
    k_point(ml::AbstractMomentumLattice, i::Int) → SVector{D, T}

The i-th k-vector in real (Cartesian) units — i.e. already multiplied
through the reciprocal basis.
"""
function k_point end

"""
    reciprocal_basis(ml::AbstractMomentumLattice) → SMatrix{D, D, T}

Reciprocal-space basis matrix. Columns are the primitive reciprocal
lattice vectors.
"""
function reciprocal_basis end

# AbstractLattice specialisations -------------------------------------
num_sites(ml::AbstractMomentumLattice) = num_k_points(ml)
position(ml::AbstractMomentumLattice, i::Int) = k_point(ml, i)

# Graph interface is vacuous for momentum lattices.
neighbors(::AbstractMomentumLattice, ::Int) = Int[]

# ---- PeriodicMomentumLattice -----------------------------------------

"""
    PeriodicMomentumLattice{D, T}

Concrete momentum lattice for a Bravais periodic lattice. Stores the
reciprocal basis, the mesh dimensions, the Γ-offset, and an explicit
list of k-points (eager construction).

Construct via [`monkhorst_pack`](@ref) or [`gamma_centered`](@ref)
rather than calling the struct directly.
"""
struct PeriodicMomentumLattice{D,T<:AbstractFloat} <: AbstractMomentumLattice{D,T}
    reciprocal_basis::SMatrix{D,D,T}
    mesh::NTuple{D,Int}
    shift::SVector{D,T}
    k_points::Vector{SVector{D,T}}
end

num_k_points(ml::PeriodicMomentumLattice) = length(ml.k_points)
k_point(ml::PeriodicMomentumLattice, i::Int) = ml.k_points[i]
reciprocal_basis(ml::PeriodicMomentumLattice) = ml.reciprocal_basis

size_trait(ml::PeriodicMomentumLattice) = FiniteSize(ml.mesh)

# Momentum lattices are implicitly periodic in k (wrapping at BZ
# boundaries). Return an all-periodic LatticeBoundary so calls to
# `boundary(ml)` don't explode.
function boundary(ml::PeriodicMomentumLattice{D}) where {D}
    axes = ntuple(_ -> PeriodicAxis(), D)
    return LatticeBoundary(axes, NoModifier())
end

periodicity(::PeriodicMomentumLattice) = Periodic()
reciprocal_support(::PeriodicMomentumLattice) = HasReciprocal()

# Momentum lattices carry no site-type information of their own. We
# give a default `UniformLayout(EmptySite())` so `site_layout` /
# `site_type` keep working if someone calls them.
site_layout(::PeriodicMomentumLattice) = UniformLayout(EmptySite())

# ---- Mesh constructors -----------------------------------------------

"""
    monkhorst_pack(basis::SMatrix{D, D, T}, mesh::NTuple{D, Int})
        → PeriodicMomentumLattice{D, T}

Construct a Monkhorst–Pack half-shifted mesh over the Brillouin zone
spanned by `basis`. For each axis the fractional k-coordinates are

    frac_i ∈ { (n + 0.5) / N - 0.5 | n = 0, …, N - 1 }

so the mesh is centred at Γ and avoids the BZ edges.
"""
function monkhorst_pack(basis::SMatrix{D,D,T}, mesh::NTuple{D,Int}) where {D,T}
    shift = zero(SVector{D,T})
    k_points = SVector{D,T}[]
    for idx in CartesianIndices(mesh)
        frac = SVector{D,T}(
            ntuple(i -> (T(idx.I[i] - 1) + T(0.5)) / T(mesh[i]) - T(0.5), D)
        )
        push!(k_points, basis * frac)
    end
    return PeriodicMomentumLattice{D,T}(basis, mesh, shift, k_points)
end

"""
    gamma_centered(basis::SMatrix{D, D, T}, mesh::NTuple{D, Int})
        → PeriodicMomentumLattice{D, T}

Construct a Γ-centred regular mesh. Fractional k-coordinates are
`n / N` with `n = 0, …, N - 1`, so the mesh includes Γ and walks to
but not through the BZ edge.
"""
function gamma_centered(basis::SMatrix{D,D,T}, mesh::NTuple{D,Int}) where {D,T}
    shift = zero(SVector{D,T})
    k_points = SVector{D,T}[]
    for idx in CartesianIndices(mesh)
        frac = SVector{D,T}(ntuple(i -> T(idx.I[i] - 1) / T(mesh[i]), D))
        push!(k_points, basis * frac)
    end
    return PeriodicMomentumLattice{D,T}(basis, mesh, shift, k_points)
end

# ---- Entry points ----------------------------------------------------

"""
    reciprocal_lattice(lat::AbstractLattice) → PeriodicMomentumLattice

Construct the reciprocal lattice for a Bravais-like `lat`.

# Trait contract

This method is a **required override** for any concrete lattice
whose `reciprocal_support(lat)` returns [`HasReciprocal`](@ref).
The fallback deliberately throws a `MethodError` so missing
implementations surface immediately at the trait dispatch site
(see [`momentum_lattice`](@ref), which routes `HasReciprocal`
lattices through `reciprocal_lattice`).

Lattices with `reciprocal_support(lat) == NoReciprocal()` MUST NOT
override this method — they have no Bravais reciprocal structure
and `momentum_lattice` will refuse to call them. Lattices with
`HasFourierModule()` should implement [`fourier_module`](@ref)
instead.

The returned object is expected to satisfy the
[`AbstractMomentumLattice`](@ref) interface (`num_k_points`,
`k_point`, `reciprocal_basis`); the canonical concrete return type
is [`PeriodicMomentumLattice`](@ref).
"""
function reciprocal_lattice(lat::AbstractLattice)
    throw(MethodError(reciprocal_lattice, (lat,)))
end

"""
    fourier_module(lat::AbstractLattice) → AbstractMomentumLattice

Quasicrystal-side entry point: construct the discrete Fourier module
(Bragg peak set) for a cut-and-project lattice.

# Trait contract

This method is a **required override** for any concrete lattice
whose `reciprocal_support(lat)` returns [`HasFourierModule`](@ref).
Concrete quasicrystal lattices implement it in their own package
(typically `QuasiCrystal.jl`); the fallback throws a `MethodError`
so the trait/method mismatch is visible at the call site through
[`momentum_lattice`](@ref).

Lattices with [`HasReciprocal`](@ref) should implement
[`reciprocal_lattice`](@ref) instead, and lattices with
[`NoReciprocal`](@ref) should not override either.

The returned object is expected to be an
[`AbstractMomentumLattice`](@ref) — typically a
[`BraggPeakSet`](@ref) — so observers like `structure_factor`
treat periodic and quasiperiodic lattices uniformly.
"""
function fourier_module(lat::AbstractLattice)
    throw(MethodError(fourier_module, (lat,)))
end

"""
    momentum_lattice(lat::AbstractLattice) → AbstractMomentumLattice

Unified trait-dispatched entry point. Returns the reciprocal lattice
for Bravais-like structures, the Fourier module for quasicrystals,
and throws for lattices without any k-space representation.
"""
momentum_lattice(lat::AbstractLattice) = _momentum_lattice(lat, reciprocal_support(lat))
_momentum_lattice(lat, ::HasReciprocal) = reciprocal_lattice(lat)
_momentum_lattice(lat, ::HasFourierModule) = fourier_module(lat)
function _momentum_lattice(lat, ::NoReciprocal)
    throw(ArgumentError("$(typeof(lat)) has no k-space representation"))
end

# ---- Quasicrystal hooks (type skeletons only) ------------------------
#
# The concrete implementations (generation algorithms, window Fourier
# transforms, symmetry reduction) live in QuasiCrystal.jl. We only
# expose the type vocabulary here so downstream code can refer to
# them.

"""
    AcceptanceWindow

Abstract supertype for acceptance windows in the internal (perp)
space of a cut-and-project quasicrystal. Concrete windows
(`PolygonalWindow`, `IntervalWindow`, etc.) live in QuasiCrystal.jl.
"""
abstract type AcceptanceWindow end

"""
    HyperReciprocalLattice{DPhys, DHyper, T}

Infinite-abstract representation of a quasicrystal's reciprocal
structure: the higher-dimensional reciprocal basis, the parallel and
perpendicular projections, and the acceptance window used to decide
peak intensities.

Concrete construction (building the projections, sampling peaks)
lives in QuasiCrystal.jl.
"""
struct HyperReciprocalLattice{DPhys,DHyper,DPerp,T<:AbstractFloat,W<:AcceptanceWindow}
    hyper_basis::SMatrix{DHyper,DHyper,T}
    parallel_proj::SMatrix{DPhys,DHyper,T}
    perp_proj::SMatrix{DPerp,DHyper,T}
    window::W
end

"""
    BraggPeakSet{DPhys, T}

Finite materialisation of a `HyperReciprocalLattice` up to a cutoff
radius: a list of physical-space k-positions with intensities and a
back-pointer to the higher-dimensional indices that generated them.

Being a subtype of [`AbstractMomentumLattice`](@ref), a `BraggPeakSet`
can be passed straight to `structure_factor` or to a
`StructureFactorObserver` so the same code paths work for periodic
and quasiperiodic lattices.
"""
struct BraggPeakSet{DPhys,DHyper,T<:AbstractFloat} <: AbstractMomentumLattice{DPhys,T}
    peaks::Vector{SVector{DPhys,T}}
    intensities::Vector{T}
    hyper_indices::Vector{NTuple{DHyper,Int}}
end

num_k_points(bps::BraggPeakSet) = length(bps.peaks)
k_point(bps::BraggPeakSet, i::Int) = bps.peaks[i]
function reciprocal_basis(::BraggPeakSet{DPhys,DHyper,T}) where {DPhys,DHyper,T}
    # A quasicrystal does not have a Bravais reciprocal basis; the
    # Fourier module is dense. We return a zero matrix as a sentinel
    # so downstream code can detect the situation.
    return zero(SMatrix{DPhys,DPhys,T})
end

size_trait(bps::BraggPeakSet) = FiniteSize((length(bps.peaks),))
function boundary(::BraggPeakSet{DPhys}) where {DPhys}
    axes = ntuple(_ -> OpenAxis(), DPhys)
    return LatticeBoundary(axes, NoModifier())
end
reciprocal_support(::BraggPeakSet) = HasFourierModule()
site_layout(::BraggPeakSet) = UniformLayout(EmptySite())
