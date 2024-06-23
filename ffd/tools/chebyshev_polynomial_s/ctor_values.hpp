namespace ffd::chebyshev_polynomial_s{

  template<std::size_t n, typename field_t, ffd::phys::statistics z,
	   template<typename, std::size_t> typename array_t>
  TPoly<n, field_t, z, array_t>::TPoly(std::array<field_t, n> const& values_at_nodes,
				       std::array<Real, 2> lower_upper_limit){
    bounds_average = .5*(lower_upper_limit[0]+lower_upper_limit[1]);
    one_over_half_diff = 2./(-lower_upper_limit[0]+lower_upper_limit[1]);


    if(!cosinus_are_computed){
      cosinus_transform = compute_cosinus_transform<n>();
      cosinus_are_computed = true;
    }


    for(std::size_t j=0; j<n; j++){
      coef[j] = 0.;
      for(std::size_t k=0; k<n; k++){
	coef[j] += values_at_nodes[k]*cosinus_transform[k+j*n];
      }
    }
  }


}//namespace
