namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   ffd::phys::statistics stat_zeta = ffd::phys::fermi,
	   template<typename, std::size_t> typename array_t = std::array>

  array_t<Complex, n>
  
  matsubara_transform(std::size_t const frequency,
		      std::size_t const integration_cuts = 10000ul){
    using ffd::core_math::Pi;
    using std::sin, std::cos;

    
    Real const delta_tau = 1./integration_cuts;
    array_t<Complex, n> matsubara_coefficients;
    TPoly<n, Real, stat_zeta, array_t> P{{0., 1.}};
    for(std::size_t j=0; j<n; ++j){
      P.coef.fill(0.);
      P.coef[j] = 1.;
      Real r_int = 0.;
      Real i_int = 0.;
      Real const omega_n = (2.*frequency+1.)*Pi;
      for(std::size_t k=0; k<integration_cuts; ++k){
	Real const low_bound = k*delta_tau;
	Real const upp_bound = delta_tau + low_bound;
	r_int += ffd::gauss_kronrod::
	  Integrate_63([=](Real t){
			 return P(t)*cos(omega_n*t);
		       }, {low_bound, upp_bound});
	i_int += ffd::gauss_kronrod::
	  Integrate_63([=](Real t){
			 return P(t)*sin(omega_n*t);
		       }, {low_bound, upp_bound});
      }
      matsubara_coefficients[j] = Complex(r_int, i_int);
    }


    return matsubara_coefficients;
  }


}//namespace
