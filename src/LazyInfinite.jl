"""
    materialize(abstract; kwargs...) → AbstractLattice

Materialise an infinite / conceptually-infinite abstract lattice into
a finite [`AbstractLattice`](@ref) up to the given cutoff.

This is the entry point for working with infinite structures
(quasicrystals beyond a cutoff radius, Fibonacci / L-system
substitutions beyond a certain depth, etc.) inside LatticeCore. It
is intentionally **generic with no typed supertype**: any package can
ship its own "infinite abstract" type and implement `materialize`
without having to inherit from a LatticeCore hierarchy.

# Cutoff convention

The meaning of the cutoff keyword argument(s) is implementation-
defined. Typical choices are:

- `depth::Int` — substitution depth (Fibonacci, L-system)
- `radius::Real` — spatial radius (Penrose cut-and-project)
- `dims::NTuple{D, Int}` — explicit per-axis sizes

The returned lattice must be [`is_finite`](@ref) `== true`, i.e. its
[`size_trait`](@ref) should be a [`FiniteSize`](@ref), so it can be
fed directly into Monte Carlo algorithms guarded by
[`require_finite`](@ref).

See `dev/note/04_architecture/06_lazy_infinite/README.md` for the
infinite-abstract ↔ finite-materialisation design pattern, and for
the parallel pattern in 05 Part B (`HyperReciprocalLattice →
BraggPeakSet`).
"""
function materialize end

"""
    require_finite(lat::AbstractLattice)

Assert that `lat` is a finite lattice. Throws `ArgumentError` if
[`is_finite`](@ref) returns `false`.

Intended as a guard at the entry point of Monte Carlo algorithms or
any routine that cannot operate on an infinite or not-yet-materialised
lattice. Returns `nothing` on success.

# Example

```julia
function run!(rng, state, lat::AbstractLattice, model, alg; kwargs...)
    require_finite(lat)
    # ... safe to walk the full site list from here ...
end
```
"""
function require_finite(lat::AbstractLattice)
    if !is_finite(lat)
        throw(
            ArgumentError(
                "$(typeof(lat)) is not finite " *
                "(size_trait = $(size_trait(lat))); call `materialize` first",
            ),
        )
    end
    return nothing
end
