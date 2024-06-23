namespace ffd::sigmadet_hubbard{

  template<typename nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real>,
	   typename G0_t=int,
	   typename Xi_t=int>
  
  auto
  compute_rho(G0_t const& G0, Xi_t const& Xi){
    int const n = ffd::core_math::sqrt_int( size(G0) );
    auto rho = Xi;
    for( int j = 0; j < n*n; ++j){
      rho[j] = nilpoly_t(n);
    }
    
    
    for( int j = 0; j < n; ++j ){
      for( int k = 0; k < n; ++k ){
	if( j != k ){
	  for( int l = 0; l < n; ++l){
	    rho[j*n+k] += G0[j*n + l]*Xi[l*n+k];
	  }
	}
      }
    }


    return rho;
  }


}//namespace
