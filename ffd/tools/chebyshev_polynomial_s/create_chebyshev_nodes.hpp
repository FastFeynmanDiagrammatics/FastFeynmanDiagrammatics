namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n>
  std::array<Real, n>
  create_cosinus_chebyshev_nodes(){
    std::array<Real, n> cosinus_cheby;
    Real const pi_over_n = ffd::core_math::Pi/n;
    for(std::size_t j=0; j<n; ++j){
      cosinus_cheby[j] = std::cos(pi_over_n*(j+0.5));
    }
    return cosinus_cheby;
  }


  template<std::size_t n>
  
  std::array<Real, n>
  
  create_chebyshev_nodes(std::array<Real, 2> const x){
    std::array<Real, n> nodes;
    auto const mean = .5*(x[0]+x[1]);
    auto const half_diff = .5*(-x[0]+x[1]);;
    auto const cos_a = create_cosinus_chebyshev_nodes<n>();
    for(std::size_t j=0; j<n; ++j){
      nodes[j] = mean+half_diff*cos_a[j];
    }
    return nodes;
  }


}//namespace
