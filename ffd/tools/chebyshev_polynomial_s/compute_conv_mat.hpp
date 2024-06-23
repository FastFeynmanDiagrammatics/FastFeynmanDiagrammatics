namespace ffd::chebyshev_polynomial_s{
  
  template<std::size_t n,
	   ffd::phys::statistics stat_zeta>
  
  std::array<Real, n*n*n>

  compute_conv_mat(){
    using poly_t = TPoly<n, Real, stat_zeta>;
    poly_t P1{{0., 1.}}, P2{{0., 1.}};
    std::array<Real, n*n*n> conv_mat;
    
    
    for(std::size_t j=0; j<n; ++j){
      P1.coef.fill(0.);
      P1.coef[j] = 1.;
      for(std::size_t m=0; m<n; ++m){
	P2.coef.fill(0.);
	P2.coef[m] = 1.;

	auto const P3 =
	  poly_t([=](Real const t){
		   return ffd::gauss_kronrod::
		     Integrate_63([=](Real t1){return P1(t1)*P2(t-t1);},
				  {0., t}) +
		     ffd::gauss_kronrod::
		     Integrate_63([=](Real t1){
				    if constexpr(stat_zeta == ffd::phys::bose){
					return P1(t1)*P2(1.+t-t1);
				      }else{
				      return -P1(t1)*P2(1.+t-t1);
				    }},
		       {t, 1.});
		 }, {0., 1.});
	for(std::size_t k=0; k<n; ++k){
	  conv_mat[m+n*(j+n*k)] = P3.coef[k];
	}
      }
    }


    return conv_mat;
  }

}//namespace
