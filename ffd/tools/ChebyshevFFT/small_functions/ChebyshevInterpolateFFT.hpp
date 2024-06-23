

namespace ffd::chebyshev_fft{

  template<typename function_t>
  std::vector<Real>
  ChebyshevInterpolateFFT(function_t F,
			  std::array<Real, 2> low_upp_bound,
			  int log2_order){
    unsigned long const pow_2_n = 1ul << log2_order;
    std::vector<Real> nodes = ReturnChebyshevNodes(low_upp_bound, log2_order);

    

    std::vector<Complex> F_values( 2*pow_2_n );
    for( std::size_t j=0; j < pow_2_n; ++j ){
      F_values[j]  =  F( nodes[j] )  ;
      F_values[j + pow_2_n]  =  0;
    }

    
    std::vector<Complex> complex_coef = ffd::fft::FFT(F_values);

    
    std::vector<Real> coef( pow_2_n );
    Real two_over_pow_2_n = 2./pow_2_n;
    using ffd::core_math::Pi;
    Real const pi_over_pow_2_np1 = -.5*Pi/pow_2_n;
    for( std::size_t k=0; k < pow_2_n; ++k ){
      coef[k] = two_over_pow_2_n *
	std::real( complex_coef[k]*
		   ffd::user_space::expI( pi_over_pow_2_np1*k ) );
    }


    return coef;
  }


}//namespace
