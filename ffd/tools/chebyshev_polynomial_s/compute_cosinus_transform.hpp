namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n>

  std::array<Real, n*n>
  compute_cosinus_transform(){
    std::array<Real, n*n> ret;


    Real constexpr pi_over_n = ffd::core_math::Pi/n;
    Real constexpr two_over_n = 2./n;
    for(std::size_t j=0; j<n; j++){
      Real const j_pi_over_n = j*pi_over_n;
      for(std::size_t k=0; k<n; k++){
	ret[k+j*n] = std::cos((k+0.5)*j_pi_over_n)*two_over_n;
      }
    }


    return ret;
  }
  


}//namespace
