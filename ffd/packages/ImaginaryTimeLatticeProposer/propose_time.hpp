#pragma once

namespace ffd::imaginary_time_lattice_proposer{

  template<typename imaginary_time_t,
	   typename lattice_t>
  Real
  Proposer<imaginary_time_t, lattice_t>::
  propose_time() const{
    Real p = ffd::user_space::Proba();

    
    auto F = [p, this](Real tau){
	       return p - time_cumulative(tau);
	     };


    Real const beta = Beta(imaginary_time_v);
    auto time = ffd::find_root::FindRealRoot(F, {-beta/2, beta/2});
    assert( time.has_value() );


    return time.value();
  }


}//namespace
