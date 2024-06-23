namespace ffd::lattice{

  template<typename lattice_t>
  
  auto
  Range(lattice_t lattice_v){
    
    
    constexpr int d = lattice_t::dimension;
    constexpr int n_atoms = lattice_t::number_atoms_unit_cell;
    constexpr int d_eff = d + (n_atoms>1);
    
    
    std::vector<decltype(CreateCoordinates(lattice_v))> vector_lattice_coordinates;
    std::vector<std::array<int, d_eff>> vector_lattice_coordinates_int;
    
    
    auto L = GetLinearSizes(lattice_v);
    
    
    for( auto r: ffd::user_space::VectorRange<d>(L) ){
      std::array<int, d_eff> components;
      for( auto j: ffd::vector_range::Range(d) ){
	components[j] = r[j];
      }
      for( [[maybe_unused]] auto atom: ffd::vector_range::Range(n_atoms) ){
	if constexpr(n_atoms > 1){
	    components[d] = atom;
	  }
	auto X = SetSpaceCoordinates(lattice_v,
				     components);
	vector_lattice_coordinates.push_back(X);
	vector_lattice_coordinates_int.push_back(components);
      }
    }

    
    return std::pair<decltype(vector_lattice_coordinates),
		     decltype(vector_lattice_coordinates_int)>{vector_lattice_coordinates,
	vector_lattice_coordinates_int};
  }
  
}//namespace
