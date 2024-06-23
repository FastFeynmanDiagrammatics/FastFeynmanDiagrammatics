namespace ffd::lattice{

  template<int coord,
	   typename lattice_t,
	   typename coordinates_t,
	   typename vector_t>
  coordinates_t
  SetSpaceCoordinatesRec(coordinates_t coordinates_v,
			 vector_t vector_v){
    ffd::get<coord>(coordinates_v).Variable = vector_v[coord];
    if constexpr( coord < lattice_t::dimension - 1 + (lattice_t::number_atoms_unit_cell > 1) ){
	return SetSpaceCoordinatesRec<coord+1, lattice_t>(coordinates_v, vector_v);
      }else{
      return coordinates_v;
    }
  }



  template<typename lattice_t,
	   typename vector_t>
  auto
  SetSpaceCoordinates(lattice_t lattice_v,
		      vector_t vector_v){
    auto coordinates_v = CreateCoordinates(lattice_v);
    return SetSpaceCoordinatesRec<0, lattice_t>(coordinates_v, vector_v);
  }
    

}//namespace
