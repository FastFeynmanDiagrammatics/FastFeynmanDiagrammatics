#pragma once

namespace ffd::imaginary_time_lattice_proposer{

  template<typename imaginary_time_t,
	   typename lattice_t>

  typename Proposer<imaginary_time_t, lattice_t>::array_space_components_t
  
  Proposer<imaginary_time_t, lattice_t>::
  propose_space(int atom) const{


    Real const p = ffd::user_space::Proba();

    
    
    auto site_iterator = std::lower_bound( begin(space_cumulative[atom]),
					   end(space_cumulative[atom]),
					   p);
    auto index_site = std::distance( begin(space_cumulative[atom]),
				     site_iterator );
    

    
    return vector_lattice_sites_components[ index_site ];

    
  }

}//namespace
