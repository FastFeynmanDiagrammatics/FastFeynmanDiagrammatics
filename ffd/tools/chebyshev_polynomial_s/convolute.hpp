namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   typename field_t,
	   ffd::phys::statistics stat_zeta,
	   template<typename, std::size_t> typename array_t>

  TPoly<n, field_t, stat_zeta, array_t>

  convolute(TPoly<n, field_t, stat_zeta, array_t> const& f,
	    TPoly<n, field_t, stat_zeta, array_t> const& g){
    using poly_t = TPoly<n, field_t, stat_zeta, array_t>;
    using conv_mat_t = TPoly_conv_mat<n, stat_zeta>;
    
    
    if( !conv_mat_t::conv_mat_is_computed ){
      conv_mat_t::conv_mat = compute_conv_mat<n, stat_zeta>();
      conv_mat_t::conv_mat_is_computed = true;
    }
    
    
    Real const factor_normalization = 2./f.one_over_half_diff;
    poly_t ret{{0., 1.}};
    for(std::size_t k=0; k<n; ++k){
      field_t coef_k = 0.;
      std::size_t const ind_k = n*n*k;
      for(std::size_t j=0; j<n; ++j){
	std::size_t const ind_shift = ind_k+n*j;
	field_t coef_k_j = 0.;
	for(std::size_t m=0; m<n; ++m){
	  coef_k_j += conv_mat_t::conv_mat[m+ind_shift]*g.coef[m];
	}
	coef_k += coef_k_j*f.coef[j];
      }
      ret.coef[k] = coef_k*factor_normalization;
    }

    
    ret.bounds_average = f.bounds_average;
    ret.one_over_half_diff = f.one_over_half_diff;
    return ret;
  }
  


}//namespace
