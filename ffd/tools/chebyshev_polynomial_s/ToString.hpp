namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   typename field_t,
	   ffd::phys::statistics stat_zeta,
	   template<typename, std::size_t> typename array_t>
  auto ToString(TPoly<n, field_t, stat_zeta, array_t> const& P,
		uint precision = 10){
    std::stringstream ss;
    ss << "ffd:TPoly" << n << ' ' << std::setprecision(precision);
    for(uint j=0; j<n; ++j){
      ss << P.coef[j] << ' ';
    }
    ss << P.bounds_average << ' ' << P.one_over_half_diff << ' ';
    return ss.str();
  }

}
