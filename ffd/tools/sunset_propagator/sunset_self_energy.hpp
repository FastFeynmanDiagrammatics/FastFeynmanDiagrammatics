namespace ffd::sunset_propagator{

  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  
  fermi_G_t<l_x, l_y, n_t>
  
  sunset_self_energy(fermi_G_t<l_x, l_y, n_t> const& G,
		     Real const U){
    fermi_G_t<l_x, l_y, n_t> Sigma;
    Real const Beta = 2.*G[0].bounds_average;
    Real const U2 = U*U;
    for(std::size_t r=0; r<size(Sigma); ++r){
      auto const G_r = G[r];
      Sigma[r] =
	fermi_cheby_s<n_t>([=](Real t){
			     Real const G_r_tau = G_r(t);
			     return U2*G_r_tau*G_r_tau*G_r(Beta-t);
			       }, ffd::make_array(0., Beta));
    }
    
    
    return Sigma;
  }

}//namespace
