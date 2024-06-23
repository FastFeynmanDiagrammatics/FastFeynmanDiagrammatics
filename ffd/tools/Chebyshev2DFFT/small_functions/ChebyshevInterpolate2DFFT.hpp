

namespace ffd::chebyshev_2d_fft{

  template<typename function2d_t>
  
  std::vector<Real>

  ChebyshevInterpolate2DFFT(function2d_t F,
			    std::array<std::array<Real, 2>, 2> low_upp_bounds,
			    std::array<unsigned long, 2> log2_sizes){
    using ulong_t = Chebyshev2DFFT::ulong_t;
    ulong_t pow_2_n = 1ul<<log2_sizes[0];
    ulong_t pow_2_m = 1ul<<log2_sizes[1];
    std::array<std::vector<Real>, 2> nodes;
    for( int j: {0, 1} ){
      nodes[j] = ffd::chebyshev_fft::
	ReturnChebyshevNodes(low_upp_bounds[j], log2_sizes[j]);
    }


    std::vector<Complex> F_values( 4*pow_2_n*pow_2_m );
    for( ulong_t y=0; y < pow_2_m; ++y ){
      for( ulong_t x=0; x < pow_2_n; ++x ){
	F_values[ x + 2*pow_2_n*y]  =  F( nodes[0][x], nodes[1][y] )  ;
	F_values[ x + pow_2_n + 2*pow_2_n*y] = 0;
	F_values[ x + pow_2_n + 2*pow_2_n*(y + pow_2_m)] = 0;
	F_values[ x + 2*pow_2_n*(y + pow_2_m)] = 0;
      }
    }

    
    std::vector<Complex> complex_coef = ffd::fft::
      FFT2D(F_values, 1ul<<(log2_sizes[0]+1) );
    std::vector<Real> coef( pow_2_n*pow_2_m );
    Real const two_over_pow_2_npm = 2./pow_2_n/pow_2_m;
    using ffd::core_math::Pi;
    Real const pi_over_pow_2_np1 = -.5*Pi/pow_2_n;
    Real const pi_over_pow_2_mp1 = -.5*Pi/pow_2_m;
    for( ulong_t kx=0; kx < pow_2_n; ++kx ){
      coef[kx] = 2.*two_over_pow_2_npm *
	std::real( complex_coef[kx]*
		   ffd::user_space::expI( pi_over_pow_2_np1 * kx ) );
    }
    for( ulong_t ky=1; ky < pow_2_m; ++ky ){
      for( ulong_t kx=0; kx < pow_2_n; ++kx ){
	ulong_t const two_2m_ky = 2*pow_2_m - ky;
	coef[kx + pow_2_n*ky] = two_over_pow_2_npm *
	  std::real( complex_coef[kx + 2*pow_2_n*ky]*
		     ffd::user_space::expI( pi_over_pow_2_np1 * kx +
					    pi_over_pow_2_mp1 * ky)
		     
		     -
		     
		     complex_coef[kx + 2*pow_2_n*two_2m_ky]*
		     ffd::user_space::expI( pi_over_pow_2_np1 * kx +
					    pi_over_pow_2_mp1 * two_2m_ky)
		     );
      }
    }

    
    return coef;
  }


}//namespace
