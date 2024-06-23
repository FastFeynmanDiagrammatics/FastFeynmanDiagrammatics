#pragma once

namespace ffd::imaginary_time_lattice_proposer{

  template<typename imaginary_time_t,
	   typename lattice_t>
  
  template<typename time_seed_function_t>
  
  void
  Proposer<imaginary_time_t, lattice_t>::
  create_time_cumulative(time_seed_function_t time_seed){


    auto const beta = Beta(imaginary_time_v);

    
    time_function = ffd::chebyshev_fft::ChebyshevFFT(time_seed,
						     {-beta/2, beta/2},
						     1e-12);

    
    
    time_cumulative = time_function.Integral();
    Real const norm = time_cumulative(beta/2);
    for(  auto  j:  ffd::vector_range::Range( size(time_function) )  ){
      time_function.Coef[j] /= norm;
      time_cumulative.Coef[j] /= norm;
    }

    
  }

}//namespace
