namespace ffd::imaginary_time_lattice{

  template<int coord,
	   typename lattice_t,
	   typename coordinates_t,
	   typename vector_t>
  
  coordinates_t
  SetSpaceCoordinatesRec(coordinates_t coordinates_v,
			 vector_t vector_v){
    ffd::get<coord+1>(coordinates_v).Variable = vector_v[coord];
    if constexpr( coord < lattice_t::dimension - 1 + (lattice_t::number_atoms_unit_cell > 1) ){
	return SetSpaceCoordinatesRec<coord+1, lattice_t>(coordinates_v, vector_v);
      }else{
      return coordinates_v;
    }
  }


  

  template<typename imaginary_time_t,
	   typename lattice_t,
	   typename pair_real_vector_t>
  
  auto
  SetCoordinates(imaginary_time_t imaginary_time_v,
		 lattice_t lattice_v,
		 pair_real_vector_t vector_v){
    auto X = ffd::user_space::
      CreateCoordinates(imaginary_time_v, lattice_v);

    
    auto [tau, x] = vector_v;
    ffd::get<0>(X).Variable = tau;


    return SetSpaceCoordinatesRec<0, lattice_t>(X, x);
  }
  


}//namespace
