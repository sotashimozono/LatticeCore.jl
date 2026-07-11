# Localization diagnostics: finite-size scaling of the inverse participation
# ratio, the standard numerical probe of Anderson / quasiperiodic
# localization. Generic over `AbstractLattice`.

"""
    mean_inverse_participation_ratio(H::AbstractMatrix; energy_window=nothing)
        -> Float64

Mean [`inverse_participation_ratio`](@ref) over the eigenstates of the real
symmetric Hamiltonian `H`. With `energy_window = (Emin, Emax)` only eigenstates
in that closed interval are averaged (mobility-edge probe); an empty window is
an error.
"""
function mean_inverse_participation_ratio(H::AbstractMatrix; energy_window=nothing)
    F = eigen(Symmetric(Matrix(H)))
    sel = _window_indices(F.values, energy_window)
    s = 0.0
    @inbounds for n in sel
        s += inverse_participation_ratio(@view F.vectors[:, n])
    end
    return s / length(sel)
end

"""
    mean_inverse_participation_ratio(lat::AbstractLattice; t=1.0, onsite=0.0,
                                     energy_window=nothing) -> Float64

Convenience method building the uniform-hopping / on-site-modulated
tight-binding Hamiltonian of `lat` and returning its mean IPR. For
bond-modulated models, build `H` with the per-bond
[`tight_binding_hamiltonian`](@ref) and pass it to the matrix form.
"""
function mean_inverse_participation_ratio(
    lat::AbstractLattice; t::Real=1.0, onsite=0.0, energy_window=nothing
)
    H = tight_binding_hamiltonian(lat; t=t, onsite=onsite)
    return mean_inverse_participation_ratio(H; energy_window=energy_window)
end

_window_indices(vals::AbstractVector, ::Nothing) = eachindex(vals)
function _window_indices(vals::AbstractVector, window)
    lo, hi = window
    sel = findall(e -> lo ≤ e ≤ hi, vals)
    isempty(sel) && throw(ArgumentError("no eigenstates in energy_window $((lo, hi))"))
    return sel
end

"""
    ipr_scaling(lats; t=1.0, onsite=0.0, energy_window=nothing)
        -> (sizes::Vector{Int}, mean_iprs::Vector{Float64})

Mean IPR of each lattice in the size-ordered collection `lats` (uniform-hopping
model), paired with `num_sites`. Feed to [`ipr_scaling_exponent`](@ref) for the
localization exponent `τ`.
"""
function ipr_scaling(lats; t::Real=1.0, onsite=0.0, energy_window=nothing)
    sizes = Int[num_sites(l) for l in lats]
    mis = Float64[
        mean_inverse_participation_ratio(
            l; t=t, onsite=onsite, energy_window=energy_window
        ) for l in lats
    ]
    return sizes, mis
end

"""
    ipr_scaling_exponent(sizes, mean_iprs) -> Float64

Localization exponent `τ` from a least-squares fit of
`log(mean_iprs) = a - τ · log(sizes)`: `τ ≈ 1` (extended), `τ ≈ 0` (localized),
`0 < τ < 1` (critical / multifractal). Needs ≥ 2 positive points.
"""
function ipr_scaling_exponent(sizes, mean_iprs)
    length(sizes) == length(mean_iprs) ||
        throw(DimensionMismatch("sizes and mean_iprs have different lengths"))
    n = length(sizes)
    n ≥ 2 || throw(ArgumentError("need at least two size points, got $n"))
    all(>(0), mean_iprs) ||
        throw(ArgumentError("all mean_iprs must be positive to take a log"))
    x = log.(float.(collect(sizes)))
    y = log.(float.(collect(mean_iprs)))
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    sxx = sum((xi - x̄)^2 for xi in x)
    sxy = sum((x[i] - x̄) * (y[i] - ȳ) for i in 1:n)
    return -(sxy / sxx)
end
