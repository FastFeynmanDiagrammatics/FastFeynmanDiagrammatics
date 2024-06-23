namespace ffd::bravais_propagator{

  template<int d,
	   typename lattice_t,
	   typename coord_t=int,
	   typename vector_t=int>
  int
  index_difference_coordinates(coord_t const& x1,
			       coord_t const& x2,
			       vector_t const& L
			       ){
    auto s1 = ffd::imaginary_time_lattice::
      GetSpaceCoordinates<lattice_t>(x1);
    auto s2 = ffd::imaginary_time_lattice::
      GetSpaceCoordinates<lattice_t>(x2);
    int index = 0, L_prod = 1;
    for(int j=0; j<d; ++j){
      int diff = s1[j]-s2[j];
      if(diff >= L[j])
	diff -= L[j];
      if(diff < 0)
	diff += L[j];
      index += L_prod*diff;
      L_prod *= L[j];
    }
    return index;
  }

}//namespace
