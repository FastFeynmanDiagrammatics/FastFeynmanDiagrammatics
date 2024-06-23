namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n, typename field_t, ffd::phys::statistics z,
	   template<typename, std::size_t> typename array_t>
  [[nodiscard]] field_t
  TPoly<n, field_t, z, array_t>::operator()(Real const x) const noexcept{
    if constexpr(n == 0){
      return 0.;
    }
    Real const t = (x - bounds_average)*one_over_half_diff;
    std::array<field_t, n+2> d_Cheby;
    d_Cheby.fill(0.);
    for(int j=int(n)-1; j >= 1; j--){
      d_Cheby[j] = d_Cheby[j+1]*2.*t - d_Cheby[j+2] + coef[j];
    }
    return t*d_Cheby[1] - d_Cheby[2] + .5*coef[0];
  }

}//namespace
