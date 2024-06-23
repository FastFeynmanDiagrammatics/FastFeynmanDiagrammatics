namespace ffd::sunset_propagator{

  template<std::size_t l_x, std::size_t l_y, std::size_t n_t,
	   std::size_t n_iter,
	   ffd::phys::diagrammatic_scheme scheme>
  
  std::array<fermi_G_t<l_x, l_y, n_t>, 2>

  third_order_propagator(fermi_G_t<l_x, l_y, n_t> const& G0,
			 Real const U,
			 Real const alpha_iter,
			 std::array<Real, 2> resummation_factors = {1., 1.}){
    static_assert(scheme == ffd::phys::bold ||
		  scheme == ffd::phys::bare ||
		  scheme == ffd::phys::renormalized_self_energy);
    auto G = G0;
    auto Sigma = sunset_self_energy<l_x, l_y, n_t>(G, U);
    auto Sigma3 = third_order_self_energy<l_x, l_y, n_t>(G, U);
    for(std::size_t r=0; r<size(Sigma); ++r){
      for(std::size_t j=0; j<n_t; ++j){
	Sigma[r].coef[j] = resummation_factors[0]*Sigma[r].coef[j]+
	  resummation_factors[1]*Sigma3[r].coef[j];
      }
    }
    
    
    if constexpr( scheme == ffd::phys::renormalized_self_energy ){
	constexpr int omega0 = 0;
	Sigma = project_to_frequency<omega0, l_x, l_y, n_t>(Sigma);
      }
    
    
    for(std::size_t iter=0; iter<n_iter; ++iter){
      G = iteration_dyson_equation<l_x, l_y, n_t>(G0, Sigma, G, alpha_iter);
      if constexpr(scheme == ffd::phys::bold ||
		   scheme == ffd::phys::renormalized_self_energy){
	  Sigma = sunset_self_energy<l_x, l_y, n_t>(G, U);
	  Sigma3 = third_order_self_energy<l_x, l_y, n_t>(G, U);
	  for(std::size_t r=0; r<size(Sigma); ++r){
	    for(std::size_t j=0; j<n_t; ++j){
	      Sigma[r].coef[j] = resummation_factors[0]*Sigma[r].coef[j]+
		resummation_factors[1]*Sigma3[r].coef[j];
	    }
	  }
	  if constexpr( scheme == ffd::phys::renormalized_self_energy ){
	      constexpr int omega0 = 0;
	      Sigma = project_to_frequency<omega0, l_x, l_y, n_t>(Sigma);
	    }
	}
    }
    
    
    return ffd::make_array(G, Sigma);
  }


}//namespace
