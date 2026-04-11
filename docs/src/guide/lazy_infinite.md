# Lazy / infinite lattices

LatticeCore has two small hooks for working with lattices that are
conceptually infinite, or that need to be materialised at a chosen
cutoff before use:

- [`materialize`](@ref) — an infinite-abstract-to-finite-lattice
  generic function.
- [`require_finite`](@ref) — a guard that throws if a lattice is
  not finite.

Together with the size trait hierarchy (see
[Traits reference](../reference/traits.md)), they cover the
lattice-side half of chapter 06 of the architecture notes. The
concept-side story lives in
[Quasiperiodic order](../concepts/quasiperiodic_order.md).

## Size traits

Every `AbstractLattice` advertises its extent through
[`size_trait(lat)`](@ref), which returns one of:

- [`FiniteSize{D}`](@ref) — an ordinary finite lattice with
  known per-axis cell counts. Safe for Monte Carlo.
- [`InfiniteSize`](@ref) — a lattice that is infinite in the
  honest sense. You can compute with it spectrally (tight binding
  on a translationally invariant chain, for instance), but you
  cannot run a Monte Carlo sweep on it.
- [`QuasiInfiniteSize{T}`](@ref) — infinite in principle, meant
  to be lowered to a finite lattice through
  [`materialize`](@ref) before being consumed. `T` is the
  numeric type of the cutoff parameter.

The predicate [`is_finite(lat)`](@ref) is a one-liner derived
from this trait:

```julia
is_finite(lat) = size_trait(lat) isa FiniteSize
```

Reference lattices always return a `FiniteSize`, so
`is_finite(LineLattice(5))` and
`is_finite(SimpleSquareLattice(3, 3))` are both `true`.

## Guarding Monte Carlo entry points

Monte Carlo algorithms walk every site, so running them on a
non-finite lattice is a bug. The canonical guard is
[`require_finite(lat)`](@ref), which throws `ArgumentError` if
[`is_finite`](@ref) is false:

```julia
function run!(rng, state, lat::AbstractLattice, model, alg; kwargs...)
    require_finite(lat)
    # ... safe to iterate 1:num_sites(lat) from here ...
end
```

The error message mentions the offending lattice's `size_trait`
and hints at `materialize`, so a user who passes an infinite
abstract lattice by mistake gets a pointer to the fix.

## The infinite-abstract → finite-materialisation pattern

The pattern LatticeCore encourages for any "infinite" structure is
to keep the abstract description as its own type and implement a
method of [`materialize`](@ref) that turns it into a `FiniteSize`
lattice. `materialize` is a generic function with **no typed
supertype**, so packages can plug in without inheriting from a
LatticeCore hierarchy.

A minimal example for a hypothetical infinite Fibonacci chain:

```julia
struct InfiniteFibonacci
    rules::Dict{Char, String}
    axiom::Char
end

function LatticeCore.materialize(inf::InfiniteFibonacci; depth::Int)
    word = inf.axiom |> Ref |> string    # start from the axiom
    for _ in 1:depth
        word = join(get(inf.rules, c, string(c)) for c in word)
    end
    n = length(word)                     # number of physical sites
    return LineLattice(n, PeriodicAxis())   # returns a FiniteSize lattice
end
```

The user's code then looks like

```julia
inf = InfiniteFibonacci(Dict('L' => "LS", 'S' => "L"), 'L')
lat = materialize(inf; depth = 10)         # LineLattice{Float64, ...}
require_finite(lat)                         # safe

run!(Random.default_rng(), initial_state(lat), lat, model, alg)
```

The same pattern applies in reciprocal space:
[`HyperReciprocalLattice`](@ref) is the "infinite abstract"
description of a cut-and-project quasicrystal's Fourier module,
and [`BraggPeakSet`](@ref) is its finite materialisation at a
cutoff. Both patterns are two halves of the same design rule:
infinite structures live as abstract types; finite slices of them
are subtypes of the main lattice hierarchy.

## Cutoff conventions

The cutoff keyword argument is **implementation-defined**. You
can make it anything that parameterises your structure, as long
as you document it clearly. Common choices:

- `depth::Int` — substitution depth (Fibonacci, L-system).
- `radius::Real` — spatial radius (Penrose cut-and-project).
- `dims::NTuple{D, Int}` — explicit per-axis sizes (generic
  rectangular sampling).

`materialize` is just a stable name; the semantics belong to each
concrete infinite abstract type.

## Unhandled cases

`materialize` has **no default implementation**, so calling it on
a type that has not specialised it raises `MethodError`:

```julia
materialize("not an infinite lattice"; depth = 1)   # MethodError
```

That is a deliberate choice: a default fallback would either
silently succeed (dangerous) or always fail with a misleading
error. `MethodError` tells the user exactly which type needs a
new method.

## See also

- [Concepts: Quasiperiodic order](../concepts/quasiperiodic_order.md)
- [Reference: Lazy / infinite](../reference/lazy_infinite.md)
- [Reference: Traits](../reference/traits.md) — size traits
