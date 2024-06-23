namespace ffd::imaginary_time_lattice_proposer{
  
  template<typename imaginary_time_t,
	   typename lattice_t>
  template<typename space_seed_function_t>
  
  void
  
  Proposer<imaginary_time_t, lattice_t>::
  create_space_cumulative(space_seed_function_t space_seed){


    constexpr int n_atoms = lattice_t::number_atoms_unit_cell;
    
    
    auto vector_lattice_sites_range = ffd::lattice::Range(lattice_v);
    auto const vector_lattice_coordinates = vector_lattice_sites_range.first;
    vector_lattice_sites_components = vector_lattice_sites_range.second;

    
    using ffd::vector_range::Range;
    for( auto atom: Range(n_atoms) ){
      auto X_atom = ffd::user_space::CreateCoordinates(imaginary_time_v, lattice_v);
      if constexpr(n_atoms > 1){
	  ffd::user_space::component<lattice_t::dimension+1>(X_atom) = atom;
	}
      auto& func  = space_function[atom];
      auto& cumul = space_cumulative[atom];
      func.resize( size(vector_lattice_sites_components) );
      cumul.resize( size(vector_lattice_sites_components) );

      
      std::size_t counter = 0;
      for( auto r: vector_lattice_coordinates ){
	auto R = RealSpace(r);
	auto XR = RealSpace(X_atom);
	for( auto j: Range(lattice_t::dimension) ){
	  R[j] -= XR[j];
	}
	func[counter] = space_seed(R);

	
	if(counter > 0){
	  cumul[counter] = cumul[counter-1];
	}else{
	  cumul[0] = 0;
	}
	cumul[counter] += func[counter];
	++counter;
      }
    

      Real const norm_space = cumul[counter-1];
      for( auto j: ffd::vector_range::Range(counter) ){
	func[j]   /= norm_space;
	cumul[j]  /= norm_space;
      }
      
    }
    
  }

}//namespace
