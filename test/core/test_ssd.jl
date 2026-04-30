using LatticeCore
using Test

# Reference SSD envelope: f(i) = sin^2(π (i + 1/2) / L) on 0-indexed
# `i ∈ 0:L-1`. Used as the ground-truth comparison for the public
# `bond_weight(::SSD, lat, i, j)` method.
ssd_ref(i, L) = sin(pi * (i + 0.5) / L)^2

# Mock lattice with InfiniteSize for the size-trait guard test. We
# define it at module scope (not inside @testset) so Aqua's piracy
# / unbound-args checks see a proper top-level type.
struct _SSDInfiniteMockLat <: LatticeCore.AbstractLattice{1,Float64} end
LatticeCore.size_trait(::_SSDInfiniteMockLat) = InfiniteSize()

@testset "SSD bond_weight: 1D LineLattice" begin
    L = 8
    lat = LineLattice(L, OpenAxis())
    ssd = SSD(float(L))

    # Axis envelope sanity: the boundary site (cx = 1, i.e. i = 0) is
    # close to zero; the half-integer center i = L/2 - 1/2 evaluates
    # to exactly 1.
    @test ssd_ref(0, L) ≈ sin(pi / 16)^2
    @test ssd_ref(0, L) ≈ 0.038060233744356625
    @test ssd_ref(L / 2 - 0.5, L) ≈ 1.0

    # Symmetry: f(i) == f(L - 1 - i).
    for i in 0:(L - 1)
        @test ssd_ref(i, L) ≈ ssd_ref(L - 1 - i, L)
    end

    # bond_weight equals the arithmetic mean of the two endpoint
    # envelopes for every nearest-neighbor bond.
    for b in bonds(lat)
        ci = b.i      # 1-based site index, also equal to the cell index
        cj = b.j
        fi = ssd_ref(ci - 1, L)
        fj = ssd_ref(cj - 1, L)
        @test bond_weight(ssd, lat, b.i, b.j) ≈ (fi + fj) / 2
    end

    # Boundary edge bond (sites 1, 2) is suppressed; the central bond
    # carries the largest weight.
    edge = bond_weight(ssd, lat, 1, 2)
    centre = bond_weight(ssd, lat, L ÷ 2, L ÷ 2 + 1)
    @test edge < 0.2
    @test centre > 0.9
    @test edge < centre

    # Monotonicity: walking from the left edge toward the centre, the
    # bond weight increases for every successive nearest-neighbor pair.
    weights = [bond_weight(ssd, lat, k, k + 1) for k in 1:(L ÷ 2 - 1)]
    @test all(weights[k] < weights[k + 1] for k in 1:(length(weights) - 1))

    # Reflection symmetry of the bond weights along the chain.
    for k in 1:(L - 1)
        @test bond_weight(ssd, lat, k, k + 1) ≈
            bond_weight(ssd, lat, L - k, L - k + 1)
    end

    # Order independence in (i, j).
    @test bond_weight(ssd, lat, 1, 2) == bond_weight(ssd, lat, 2, 1)
end

@testset "SSD bond_weight: 2D SimpleSquareLattice (axis product)" begin
    Lx, Ly = 4, 6
    lat = SimpleSquareLattice(Lx, Ly, OpenAxis())
    ssd = SSD(float(min(Lx, Ly)))

    # 2D envelope at site (cx, cy) (1-based) is the product of the
    # per-axis envelopes. Verify on every site via the self-bond
    # mean: bond_weight(SSD, lat, i, i) is the envelope itself.
    for i in 1:num_sites(lat)
        x = mod1(i, Lx)
        y = ((i - 1) ÷ Lx) + 1
        f_expected = ssd_ref(x - 1, Lx) * ssd_ref(y - 1, Ly)
        @test bond_weight(ssd, lat, i, i) ≈ f_expected
    end

    # Real bond weight: arithmetic mean of two endpoint envelopes.
    for b in bonds(lat)
        x_i = mod1(b.i, Lx)
        y_i = ((b.i - 1) ÷ Lx) + 1
        x_j = mod1(b.j, Lx)
        y_j = ((b.j - 1) ÷ Lx) + 1
        fi = ssd_ref(x_i - 1, Lx) * ssd_ref(y_i - 1, Ly)
        fj = ssd_ref(x_j - 1, Lx) * ssd_ref(y_j - 1, Ly)
        @test bond_weight(ssd, lat, b.i, b.j) ≈ (fi + fj) / 2
    end

    # Corner bonds are most suppressed; an interior bulk bond carries
    # a much larger weight.
    corner_bond_weight = bond_weight(ssd, lat, 1, 2)             # along x at corner
    bulk_i = (Ly ÷ 2 - 1) * Lx + Lx ÷ 2                          # near-centre site
    bulk_bond_weight = bond_weight(ssd, lat, bulk_i, bulk_i + 1)
    @test corner_bond_weight < bulk_bond_weight
end

@testset "SSD bond_weight: BC-agnostic interface" begin
    # The implementation evaluates the envelope regardless of the
    # axis tuple. Bulk PBC and OBC must agree on the weight value
    # because the weight depends only on lattice cell coordinates,
    # not on which bonds the BC connects.
    L = 6
    obc = LineLattice(L, OpenAxis())
    pbc = LineLattice(L, PeriodicAxis())
    ssd = SSD(float(L))

    for k in 1:(L - 1)
        @test bond_weight(ssd, obc, k, k + 1) ≈ bond_weight(ssd, pbc, k, k + 1)
    end
end

@testset "SSD bond_weight: rejects non-finite size" begin
    @test_throws ArgumentError bond_weight(SSD(8.0), _SSDInfiniteMockLat(), 1, 2)
end

@testset "SSD bond_weight: type stability and value type" begin
    lat = LineLattice(4, OpenAxis())
    ssd = SSD(4.0)
    w = bond_weight(ssd, lat, 1, 2)
    @test w isa Float64
    @test w >= 0.0
    @test w <= 1.0
end
