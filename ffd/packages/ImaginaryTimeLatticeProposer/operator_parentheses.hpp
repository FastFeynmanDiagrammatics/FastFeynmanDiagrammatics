namespace ffd::imaginary_time_lattice_proposer{

  template<typename imaginary_time_t,
	   typename lattice_t>

  Real

  Proposer<imaginary_time_t, lattice_t>::
  operator()(typename Proposer<imaginary_time_t, lattice_t>::coordinates_t X1,
	     typename Proposer<imaginary_time_t, lattice_t>::coordinates_t X2) const{
    using ffd::vector_range::Range;

    
    auto [tau1, r1] = ffd::imaginary_time_lattice::
      GetCoordinates<imaginary_time_t, lattice_t>(X1);
    auto [tau2, r2] = ffd::imaginary_time_lattice::
      GetCoordinates<imaginary_time_t, lattice_t>(X2);


    Real const beta = Beta(imaginary_time_v);
    Real tau1_2 = tau1 - tau2;
    while( tau1_2 < -beta/2 ){
      tau1_2 += beta;
    }
    while( tau1_2 > beta/2 ){
      tau1_2 -= beta;
    }
    
    
    auto L = ffd::lattice::GetLinearSizes( lattice_v );
    array_space_components_t r1_2 = r1;
    for(  auto  j:  Range( lattice_t::dimension )  ){
      r1_2[j] -= r2[j];
      r1_2[j]  =  ( r1_2[j] >=0 ? r1_2[j] : r1_2[j] += L[j]);
    }
    int atom = 0;
    constexpr bool has_unit_cell = (lattice_t::number_atoms_unit_cell>1);
    if constexpr(has_unit_cell){
	atom = r2[lattice_t::dimension];
      }
    

    static_assert( std::is_same_v<decltype(r1_2), array_space_components_t> );
    
    
    auto site_iterator =
      std::find(
      // std::lower_bound(
		       begin(vector_lattice_sites_components),
		       end(vector_lattice_sites_components),
		       r1_2
		       );
    auto site_index = std::distance( begin(vector_lattice_sites_components),
    				     site_iterator );
   
    
    return time_function(tau1_2) * space_function[atom][site_index];
  }

}//namespace
