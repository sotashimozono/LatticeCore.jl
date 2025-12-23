Base.size(lat::AbstractLattice) = error("Not implemented for type $(typeof(lat)))")
boundary_conditions(lat::AbstractLattice) = error("Not implemented for type $(typeof(lat)))")

site_coordinate(lat::AbstractLattice, i::Int) = error("Not implemented for type $(typeof(lat)))")
lattice_center(lat::AbstractLattice) = sum(site_coordinate(lat, i) for i in 1:length(lat)) / length(lat)

lazy_bond_coords(lat) = (bond_coordinate(lat, b) for b in bonds(lat))