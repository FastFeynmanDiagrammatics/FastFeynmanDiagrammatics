namespace ffd::imaginary_time_lattice_proposer{

  
  template<typename imaginary_time_t,
	   typename lattice_t>

  
  typename Proposer<imaginary_time_t, lattice_t>::coordinates_t

  
  Proposer<imaginary_time_t, lattice_t>::
  propose(typename Proposer<imaginary_time_t, lattice_t>::coordinates_t X) const{
    Real const beta = Beta(imaginary_time_v);
    Real delta_tau = propose_time() + ffd::get<0>(X)();
    while( delta_tau > beta){
      delta_tau -= beta;
    }
    while( delta_tau < 0 ){
      delta_tau += beta;
    }

    
    int atom = 0;
    constexpr bool has_unit_cell = (lattice_t::number_atoms_unit_cell > 1);
    if constexpr( has_unit_cell ){
	atom = ffd::get<lattice_t::dimension+1>(X)();
      }
    auto delta_x = propose_space(atom);
    auto space_x = ffd::imaginary_time_lattice::
      GetSpaceCoordinates<lattice_t>(X);
    for(  auto  j:  ffd::vector_range::Range( lattice_t::dimension )  ){
      delta_x[j] += space_x[j];
    }


    
    return ffd::imaginary_time_lattice::
      SetCoordinates(imaginary_time_v, lattice_v, std::pair<Real,
		     decltype(delta_x)>{delta_tau, delta_x});
  }


}//namespace
