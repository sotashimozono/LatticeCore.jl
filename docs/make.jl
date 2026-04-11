using LatticeCore
using Documenter
using Downloads

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)
favicon_path = joinpath(assets_dir, "favicon.ico")
logo_path = joinpath(assets_dir, "logo.png")

# Best-effort favicon / logo download. Skipped gracefully if the
# network is unavailable (e.g. local sandboxed builds).
try
    Downloads.download("https://github.com/sotashimozono.png", favicon_path)
    Downloads.download("https://github.com/sotashimozono.png", logo_path)
catch err
    @warn "Could not download favicon/logo; continuing without them" exception = err
end

makedocs(;
    sitename = "LatticeCore.jl",
    format = Documenter.HTML(;
        canonical = "https://codes.sota-shimozono.com/LatticeCore.jl/stable/",
        prettyurls = get(ENV, "CI", "false") == "true",
        mathengine = MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
                    :tags => "ams",
                    :packages => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
        assets = ["assets/favicon.ico", "assets/custom.css"],
    ),
    modules = [LatticeCore],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Concepts" => [
            "Overview" => "concepts/index.md",
            "Lattices and unit cells" => "concepts/lattice_and_unit_cells.md",
            "Reciprocal lattice and Brillouin zone" => "concepts/reciprocal_and_brillouin.md",
            "Boundary conditions" => "concepts/boundary_conditions.md",
            "Site types and spin models" => "concepts/site_types_and_spins.md",
            "Quasiperiodic order" => "concepts/quasiperiodic_order.md",
        ],
        "Guide" => [
            "Lattice interface" => "guide/lattice.md",
            "Boundary conditions" => "guide/boundary.md",
            "Coordinate systems" => "guide/coordinates.md",
            "Site types and layouts" => "guide/site_type.md",
            "Momentum space" => "guide/momentum.md",
            "Lazy / infinite lattices" => "guide/lazy_infinite.md",
        ],
        "API Reference" => [
            "Lattice" => "reference/lattice.md",
            "Bond" => "reference/bond.md",
            "Traits" => "reference/traits.md",
            "Boundary" => "reference/boundary.md",
            "Coordinates" => "reference/coordinates.md",
            "Indexing" => "reference/indexing.md",
            "Site type" => "reference/site_type.md",
            "Site layout" => "reference/site_layout.md",
            "Lattice element" => "reference/lattice_element.md",
            "Momentum space" => "reference/momentum.md",
            "Lazy / infinite" => "reference/lazy_infinite.md",
            "Reference lattices" => "reference/reference_lattices.md",
        ],
        "Design notes" => "design.md",
    ],
    checkdocs = :none,
)

deploydocs(; repo = "github.com/sotashimozono/LatticeCore.jl.git", devbranch = "main")
