namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   typename field_t = Real,
	   ffd::phys::statistics stat_zeta = ffd::phys::fermi,
	   template<typename, std::size_t> typename array_t = std::array>
  
  Real
  norm_max(TPoly<n, field_t, stat_zeta, array_t> const& P){
    Real norm = 0.;
    for(int j=0; j<n; ++j){
      norm = std::max(std::abs(P.coef[j]), norm);
    }
    return norm;
  }


}//namespace
