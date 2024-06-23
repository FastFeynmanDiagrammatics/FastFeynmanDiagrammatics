namespace ffd::sunset_propagator{

  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  
  fermi_G_t<l_x, l_y, n_t>

  iteration_dyson_equation(fermi_G_t<l_x, l_y, n_t> const& G0,
			   fermi_G_t<l_x, l_y, n_t> const& Sigma,
			   fermi_G_t<l_x, l_y, n_t> const& G,
			   Real const alpha_iter){
    auto const G0_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G0);
    auto const Sigma_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(Sigma);
    auto const G_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(G);
    std::size_t constexpr size_x_y = (1ul<<(l_x+l_y));


    fermi_G_t<l_x, l_y, n_t> G_new_kt;
    for(std::size_t k=0; k<size_x_y; ++k){
      auto const G0_Sigma = convolute(G0_kt[k], Sigma_kt[k]);
      auto const G0_Sigma_G = convolute(G0_Sigma, G_kt[k]);
      G_new_kt[k].bounds_average = G_kt[k].bounds_average;
      G_new_kt[k].one_over_half_diff = G_kt[k].one_over_half_diff;
      for(std::size_t j=0; j<n_t; ++j){
	G_new_kt[k].coef[j] = alpha_iter*G_kt[k].coef[j]+
	  (1-alpha_iter)*(G0_kt[k].coef[j] + G0_Sigma_G.coef[j]);
      }
    }


    auto G_new_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(G_new_kt);

    
    return G_new_rt;
  }

}//namespace
