# Single-particle tight-binding models and spectral tools on any
# `AbstractLattice`. These are pure functions of the declared connectivity /
# geometry (via `adjacency_matrix`, `bonds`, `position`), so every concrete
# lattice — periodic (`Lattice2D`) or aperiodic (`QuasiCrystal`) — gets them.

"""
    tight_binding_hamiltonian(lat::AbstractLattice; t=1.0, onsite=0.0)
        -> SparseMatrixCSC{Float64,Int}

Real symmetric single-particle tight-binding Hamiltonian of `lat` with uniform
hopping `t` on every nearest-neighbour bond,

```math
H = \\mathrm{diag}(ε) - t \\, A,
```

where `A` is the [`adjacency_matrix`](@ref). `onsite` is a scalar on-site
energy or a length-`num_sites(lat)` vector. Requires a finite lattice.

See also [`spectrum`](@ref), [`inverse_participation_ratio`](@ref),
[`density_of_states`](@ref).
"""
function tight_binding_hamiltonian(lat::AbstractLattice; t::Real=1.0, onsite=0.0)
    N = num_sites(lat)
    H = (-float(t)) .* adjacency_matrix(lat; sparse=true)
    ov = _onsite_vector(onsite, N)
    if any(!iszero, ov)
        H = H + spdiagm(0 => ov)
    end
    return H
end

"""
    tight_binding_hamiltonian(lat::AbstractLattice, hoppings::AbstractVector;
                              onsite=0.0) -> SparseMatrixCSC{Float64,Int}

Per-bond hopping variant: `hoppings[k]` is the amplitude on the `k`-th bond of
`collect(bonds(lat))`. Use this for bond-modulated models (e.g. a Fibonacci
`tL`/`tS` chain). `length(hoppings)` must equal the number of bonds.
"""
function tight_binding_hamiltonian(
    lat::AbstractLattice, hoppings::AbstractVector; onsite=0.0
)
    bs = collect(bonds(lat))
    length(hoppings) == length(bs) || throw(
        DimensionMismatch(
            "hoppings has length $(length(hoppings)) but lattice has $(length(bs)) bonds",
        ),
    )
    N = num_sites(lat)
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    sizehint!(Is, 2 * length(bs) + N)
    sizehint!(Js, 2 * length(bs) + N)
    sizehint!(Vs, 2 * length(bs) + N)
    @inbounds for k in eachindex(bs)
        b = bs[k]
        tij = -float(hoppings[k])
        push!(Is, b.i);
        push!(Js, b.j);
        push!(Vs, tij)
        push!(Is, b.j);
        push!(Js, b.i);
        push!(Vs, tij)
    end
    ov = _onsite_vector(onsite, N)
    @inbounds for i in 1:N
        ov[i] != 0 && (push!(Is, i); push!(Js, i); push!(Vs, ov[i]))
    end
    return sparse(Is, Js, Vs, N, N)
end

_onsite_vector(onsite::Real, N::Int) = fill(Float64(onsite), N)
function _onsite_vector(onsite::AbstractVector, N::Int)
    length(onsite) == N ||
        throw(DimensionMismatch("onsite vector has length $(length(onsite)), expected $N"))
    return Float64.(onsite)
end

"""
    spectrum(lat::AbstractLattice; t=1.0, onsite=0.0) -> Vector{Float64}

Sorted eigenvalues of the tight-binding Hamiltonian of `lat`.
"""
function spectrum(lat::AbstractLattice; t::Real=1.0, onsite=0.0)
    return eigvals(Symmetric(Matrix(tight_binding_hamiltonian(lat; t=t, onsite=onsite))))
end

"""
    eigenstates(lat::AbstractLattice; t=1.0, onsite=0.0)
        -> (values::Vector{Float64}, vectors::Matrix{Float64})

Full eigendecomposition of the tight-binding Hamiltonian: `values` ascending,
`vectors[:, n]` the eigenvector for `values[n]`.
"""
function eigenstates(lat::AbstractLattice; t::Real=1.0, onsite=0.0)
    F = eigen(Symmetric(Matrix(tight_binding_hamiltonian(lat; t=t, onsite=onsite))))
    return F.values, F.vectors
end

"""
    inverse_participation_ratio(ψ::AbstractVector) -> Float64

Inverse participation ratio `Σ|ψ_i|⁴ / (Σ|ψ_i|²)²` of a single state. Ranges in
`[1/N, 1]`: `1/N` for a fully extended state, `1` for a state on one site; the
participation number `1/IPR` estimates how many sites it occupies.
"""
function inverse_participation_ratio(ψ::AbstractVector)
    n2 = zero(real(eltype(ψ)))
    n4 = zero(real(eltype(ψ)))
    @inbounds for a in ψ
        p = abs2(a)
        n2 += p
        n4 += p * p
    end
    n2 == 0 && throw(ArgumentError("inverse_participation_ratio of the zero vector"))
    return n4 / (n2 * n2)
end

"""
    inverse_participation_ratios(lat::AbstractLattice; t=1.0, onsite=0.0)
        -> (energies::Vector{Float64}, iprs::Vector{Float64})

Per-eigenstate IPR of the tight-binding spectrum of `lat`, with energies.
"""
function inverse_participation_ratios(lat::AbstractLattice; t::Real=1.0, onsite=0.0)
    vals, vecs = eigenstates(lat; t=t, onsite=onsite)
    iprs = [inverse_participation_ratio(@view vecs[:, n]) for n in axes(vecs, 2)]
    return vals, iprs
end

"""
    density_of_states(lat::AbstractLattice; t=1.0, onsite=0.0, nbins=100,
                      broadening=0.0) -> (centers, dos)

Density of states of the tight-binding spectrum, normalized so that
`sum(dos) * step ≈ num_sites(lat)`. `broadening = 0` gives a histogram over
`nbins` bins; `broadening = σ > 0` smears each level by a Gaussian of width `σ`.
"""
function density_of_states(
    lat::AbstractLattice; t::Real=1.0, onsite=0.0, nbins::Int=100, broadening::Real=0.0
)
    nbins ≥ 1 || throw(ArgumentError("nbins must be ≥ 1, got $nbins"))
    E = spectrum(lat; t=t, onsite=onsite)
    return _density_of_states(E, nbins, Float64(broadening))
end

function _density_of_states(E::AbstractVector, nbins::Int, σ::Float64)
    Emin, Emax = extrema(E)
    if Emax ≈ Emin
        Emin -= 0.5
        Emax += 0.5
    end
    if σ > 0
        Emin -= 3σ
        Emax += 3σ
    end
    edges = range(Emin, Emax; length=nbins + 1)
    dE = step(edges)
    centers = collect(edges[1:(end - 1)] .+ dE / 2)
    dos = zeros(Float64, nbins)
    if σ > 0
        norm = 1 / (σ * sqrt(2π))
        @inbounds for e in E, b in 1:nbins
            dos[b] += norm * exp(-((centers[b] - e)^2) / (2σ^2))
        end
    else
        @inbounds for e in E
            b = clamp(floor(Int, (e - Emin) / dE) + 1, 1, nbins)
            dos[b] += 1
        end
        dos ./= dE
    end
    return centers, dos
end
