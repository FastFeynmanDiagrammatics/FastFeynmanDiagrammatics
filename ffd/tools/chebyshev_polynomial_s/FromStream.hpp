namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n,
	   typename field_t,
	   ffd::phys::statistics stat_zeta,
	   template<typename, std::size_t> typename array_t>
  template<class stream_t>
  void TPoly<n, field_t, stat_zeta, array_t>::
  FromStream(stream_t& stream){
    std::string discard;
    stream >> discard;
    for(uint j=0; j<n; ++j){
      stream >> this->coef[j];
    }
    stream >> this->bounds_average;
    stream >> this->one_over_half_diff;
  }

}
