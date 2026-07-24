# Scale changes: generating a size sequence of lattices of the same family.
#
# Trait types live in `Traits.jl`; this file holds the accessors and the generic
# fallbacks. Concrete lattice packages provide `rescale` / `cell_partition`.

"""
    scaling_rule(lat::AbstractLattice) -> AbstractScalingRule

How this lattice family changes scale. Defaults to [`NoScaling`](@ref), i.e. the
lattice offers no canonical size sequence and [`rescale`](@ref) is unavailable.

A Bravais-like lattice normally reports [`LinearScaling`](@ref) — its per-axis
cell counts simply multiply. An aperiodic lattice normally reports
[`SubstitutionScaling`](@ref), where the natural sequence of sizes is generated
by its substitution rule rather than by an integer side length (Fibonacci
lengths, Penrose inflation radii, …), so "the next size up" is a well defined
notion even though there is no `L` to double.

See also [`rescale`](@ref), [`size_sequence`](@ref), [`cell_partition`](@ref).
"""
scaling_rule(::AbstractLattice) = NoScaling()

"""
    rescale(lat::AbstractLattice, k::Integer = 1) -> AbstractLattice

The same lattice `k` scale steps larger, in the sense of [`scaling_rule`](@ref).
`k = 0` returns a lattice equivalent to `lat`; negative `k` steps down where the
family supports it.

For [`LinearScaling`](@ref)`(f)` one step multiplies every per-axis cell count by
`f`. For [`SubstitutionScaling`](@ref)`(d)` one step advances the substitution
depth by `d`. The point of the common verb is that a routine which studies how
some quantity behaves as the system grows can be written once and applied to
periodic and aperiodic lattices alike, instead of threading `(Lx, Ly)` through
code that a quasicrystal cannot supply.

Concrete lattice packages implement this; there is no meaningful generic
fallback, so the default throws.
"""
function rescale(lat::AbstractLattice, k::Integer=1)
    return throw(
        ArgumentError(
            "rescale is not implemented for $(typeof(lat)) (scaling_rule = $(scaling_rule(lat))). " *
            "A lattice family must define `rescale(lat, k)` to take part in size sequences.",
        ),
    )
end

"""
    size_sequence(lat::AbstractLattice, n::Integer) -> Vector{<:AbstractLattice}

`[rescale(lat, k) for k in 0:n]` — the lattice at successive scales, for
sweeping a quantity against system size without hard-coding how the family is
parameterized.

```julia
for l in size_sequence(lat, 4)
    push!(xs, num_sites(l))
    push!(ys, measure_something(l))
end
```
"""
size_sequence(lat::AbstractLattice, n::Integer) = [rescale(lat, k) for k in 0:n]

"""
    cell_partition(lat::AbstractLattice, k::Integer = 1) -> Vector{Vector{Int}}

Group the sites of `lat` by which cell of the lattice `k` steps **coarser** —
that is, of `rescale(lat, -k)` — they fall into. Entry `c` lists the site indices
of `lat` belonging to cell `c` of that coarser lattice, and together the groups
partition `1:num_sites(lat)`.

Note the direction: `rescale(lat, k)` goes *up* in size and has more cells than
`lat`, so it is the `-k` lattice whose cells are unions of `lat`'s. A family that
cannot take the downward step (`LinearScaling` with cell counts not divisible by
`factor^k`, say) should raise rather than round.

Useful whenever a quantity defined per site has to be compared across scales —
cell averages of a local observable, cell-resolved densities, or transferring a
configuration between two members of a [`size_sequence`](@ref).

Concrete lattice packages implement this; the default throws.
"""
function cell_partition(lat::AbstractLattice, k::Integer=1)
    return throw(
        ArgumentError(
            "cell_partition is not implemented for $(typeof(lat)) " *
            "(scaling_rule = $(scaling_rule(lat))).",
        ),
    )
end
