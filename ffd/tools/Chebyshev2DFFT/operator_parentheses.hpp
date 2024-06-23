

namespace ffd::chebyshev_2d_fft{

  Real

  Chebyshev2DFFT::
  operator()(std::array<Real, 2> X) const{

    
    std::array<Real, 2> t;
    for(int j: {0, 1}){
      auto const [half_mean, half_diff] = ffd::chebyshev_fft::
	return_half_mean_diff(LowerUpperBounds[j]);
      t[j] = (X[j] - half_mean) / half_diff;
    }


    auto const [pow_2_n, pow_2_m] = sizes(*this);

    
    std::vector<Real> f_x_ky( pow_2_m );
    for(ulong_t ky = 0; ky < pow_2_m; ++ky){
      std::vector<Real> recursion_vec(pow_2_n+2);
      recursion_vec[pow_2_n+1] = 0;
      recursion_vec[pow_2_n] = 0;
      for(ulong_t kx=pow_2_n-1; kx >= 1; --kx){
	recursion_vec[kx] = 2*t[0]*recursion_vec[kx+1] - recursion_vec[kx+2] +
	  Coef[kx+pow_2_n*ky];
      }
      f_x_ky[ky] = t[0]*recursion_vec[1] - recursion_vec[2] +
	0.5*Coef[0+pow_2_n*ky];
    }

    
    std::vector<Real> recursion_vec(pow_2_m+2);
    recursion_vec[pow_2_m+1] = 0;
    recursion_vec[pow_2_m] = 0;
    for(ulong_t ky = pow_2_m - 1; ky >= 1; --ky){
      recursion_vec[ky] = 2*t[1]*recursion_vec[ky+1] - recursion_vec[ky+2] +
	f_x_ky[ky];
    }
    return t[1]*recursion_vec[1] - recursion_vec[2] +
      0.5*f_x_ky[0];
  }

}//namespace
