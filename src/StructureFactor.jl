"""
    structure_factor(lat, state, k::SVector) → Float64

Scalar structure factor at a single k-point:

    S(k) = (1/N) |⟨ Σ_i s_i exp(-i k · r_i) ⟩|²

computed on a single snapshot (no thermal averaging). The default
implementation is naive O(N) per k-point; FFT-based specialisations
on regular meshes can be added later without changing the API.
"""
function structure_factor(
    lat::AbstractLattice, state::AbstractVector, k::SVector{D,T}
) where {D,T}
    N = num_sites(lat)
    acc = zero(ComplexF64)
    for i in 1:N
        acc += state[i] * cis(-dot(k, position(lat, i)))
    end
    return abs2(acc) / N
end

"""
    structure_factor(lat, state, ml::AbstractMomentumLattice) → Vector{Float64}

Evaluate the structure factor at every k-point of `ml`. Naive
O(N · M) where `M = num_k_points(ml)`.
"""
function structure_factor(
    lat::AbstractLattice, state::AbstractVector, ml::AbstractMomentumLattice
)
    out = Vector{Float64}(undef, num_k_points(ml))
    for i in 1:num_k_points(ml)
        out[i] = structure_factor(lat, state, k_point(ml, i))
    end
    return out
end
