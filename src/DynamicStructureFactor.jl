# Single-particle momentum-resolved spectral function A(k,ω), the dynamic
# counterpart of the static `structure_factor`. Generic over `AbstractLattice`.

"""
    dynamic_structure_factor(lat::AbstractLattice,
                             kpoints::AbstractVector{<:AbstractVector};
                             t=1.0, onsite=0.0, omega=nothing, nomega=200,
                             broadening=0.05) -> (omegas, A)

Momentum-resolved single-particle spectral function of the tight-binding model
on `lat`,

```math
A(k, ω) = \\sum_n |\\langle k | ψ_n \\rangle|^2 \\, δ(ω - E_n),
\\qquad
\\langle k | ψ_n \\rangle = \\tfrac{1}{\\sqrt N} \\sum_l e^{-i k · r_l} ψ_n(l),
```

evaluated at every k-vector in `kpoints`. Returns the shared `omegas` grid
(length `nomega`) and a `length(kpoints) × nomega` matrix `A`. Each δ-peak is a
normalized Gaussian of width `broadening`, so `sum(A[j, :]) * Δω ≈ 1` for every
`j` (spectral-weight sum rule). On a periodic crystal `A(k,ω)` collapses onto
the Bloch bands; on a quasicrystal it resolves the fragmented spectrum.

See also [`structure_factor`](@ref), [`spectrum`](@ref).
"""
function dynamic_structure_factor(
    lat::AbstractLattice,
    kpoints::AbstractVector{<:AbstractVector};
    t::Real=1.0,
    onsite=0.0,
    omega=nothing,
    nomega::Int=200,
    broadening::Real=0.05,
)
    nomega ≥ 1 || throw(ArgumentError("nomega must be ≥ 1, got $nomega"))
    broadening > 0 || throw(ArgumentError("broadening must be > 0, got $broadening"))

    E, V = eigenstates(lat; t=t, onsite=onsite)
    N = length(E)
    R = _positions_matrix(lat, N)

    omegas = _omega_grid(omega, E, Float64(broadening), nomega)
    σ = Float64(broadening)
    gnorm = 1 / (σ * sqrt(2π))

    nk = length(kpoints)
    A = zeros(Float64, nk, length(omegas))
    phase = Vector{ComplexF64}(undef, N)
    @inbounds for j in 1:nk
        k = kpoints[j]
        for l in 1:N
            kr = 0.0
            for d in eachindex(k)
                kr += k[d] * R[d, l]
            end
            phase[l] = cis(-kr)
        end
        for n in 1:N
            c = zero(ComplexF64)
            for l in 1:N
                c += phase[l] * V[l, n]
            end
            w = abs2(c) / N
            En = E[n]
            for m in eachindex(omegas)
                A[j, m] += w * gnorm * exp(-((omegas[m] - En)^2) / (2σ^2))
            end
        end
    end
    return omegas, A
end

"""
    dynamic_structure_factor(lat::AbstractLattice, ml::AbstractMomentumLattice;
                             kwargs...)

Take the k-points from a momentum lattice (e.g. a `BraggPeakSet`), matching the
grid used by the static [`structure_factor`](@ref).
"""
function dynamic_structure_factor(
    lat::AbstractLattice, ml::AbstractMomentumLattice; kwargs...
)
    kpoints = [k_point(ml, j) for j in 1:num_k_points(ml)]
    return dynamic_structure_factor(lat, kpoints; kwargs...)
end

function _omega_grid(omega, E::AbstractVector, σ::Float64, nomega::Int)
    omega === nothing || return collect(Float64.(omega))
    Emin, Emax = extrema(E)
    if Emax ≈ Emin
        Emin -= 0.5
        Emax += 0.5
    end
    return collect(range(Emin - 3σ, Emax + 3σ; length=nomega))
end

_spatial_dim(::AbstractLattice{D}) where {D} = D

function _positions_matrix(lat::AbstractLattice, N::Int)
    D = _spatial_dim(lat)
    R = Matrix{Float64}(undef, D, N)
    @inbounds for l in 1:N
        p = position(lat, l)
        for d in 1:D
            R[d, l] = Float64(p[d])
        end
    end
    return R
end
