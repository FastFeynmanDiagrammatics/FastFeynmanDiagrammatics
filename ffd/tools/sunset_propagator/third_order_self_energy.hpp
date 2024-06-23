namespace ffd::sunset_propagator{

  template<std::size_t l_x, std::size_t l_y, std::size_t n_t>
  
  fermi_G_t<l_x, l_y, n_t>
  
  third_order_self_energy(fermi_G_t<l_x, l_y, n_t> const& G,
			  Real const U){
    fermi_G_t<l_x, l_y, n_t> Sigma;
    Real const Beta = 2.*G[0].bounds_average;
    Real const mU3 = -U*U*U;
    std::size_t constexpr size_x_y = (1ul<<(l_x+l_y));
    fermi_G_t<l_x, l_y, n_t> ph_bubble_rt, pp_bubble_rt;
    
    
    for(std::size_t r=0; r<size_x_y; ++r){
      auto const G_r = G[r];
      ph_bubble_rt[r] =
	fermi_cheby_s<n_t>([=](Real t){
			     return -G_r(t)*G_r(Beta-t);
			   }, ffd::make_array(0., Beta));
      pp_bubble_rt[r] =
	fermi_cheby_s<n_t>([=](Real t){
			     return G_r(t)*G_r(t);
			   }, ffd::make_array(0., Beta));
    }
    auto ph_bubble_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(ph_bubble_rt);
    auto pp_bubble_kt = from_spacetime_to_ktime<l_x, l_y, n_t>(pp_bubble_rt);

    
    fermi_G_t<l_x, l_y, n_t> ph_bubbles_kt, pp_bubbles_kt;
    for(std::size_t k=0; k<size_x_y; ++k){
      ph_bubbles_kt[k] = convolute(ph_bubble_kt[k], ph_bubble_kt[k]);
      pp_bubbles_kt[k] = convolute(pp_bubble_kt[k], pp_bubble_kt[k]);
    }
    

    auto ph_bubbles_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(ph_bubbles_kt);
    auto pp_bubbles_rt = from_ktime_to_spacetime<l_x, l_y, n_t>(pp_bubbles_kt);


    for(std::size_t r=0; r<size_x_y; ++r){
      Sigma[r] =
	fermi_cheby_s<n_t>([=](Real t){
			     return mU3*(G[r](t)*ph_bubbles_rt[r](t)+
					-G[r](Beta-t)*pp_bubbles_rt[r](t));
			       }, ffd::make_array(0., Beta));
    }
    
    
    return Sigma;
  }


}//namespace
