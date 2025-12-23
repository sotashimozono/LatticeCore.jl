# LatticeCore.jl の中に追加するイメージ

"""return the coordinate of site i in the lattice lat."""
function site_coordinate(lat::AbstractLattice, i::Int)
    # 実装は Lattice2D や QuasiCrystal 側で行う
    error("site_coordinate not implemented for $(typeof(lat))")
end

"""
    bond_coordinate(lat::AbstractLattice, bond::Bond)
ボンド（結合）の中心座標を返す。デフォルトでは始点と終点の中点。
"""
function bond_coordinate(lat::AbstractLattice, bond::Bond)
    r_src = site_coordinate(lat, bond.src)
    r_dst = site_coordinate(lat, bond.dst)
    return (r_src .+ r_dst) ./ 2
end

"""
    bond_coordinate(lat::AbstractLattice, i::Int, j::Int)
サイト i と j を結ぶ結合の中心座標を返す。
"""
function bond_coordinate(lat::AbstractLattice, i::Int, j::Int)
    return (site_coordinate(lat, i) .+ site_coordinate(lat, j)) ./ 2
end