namespace ffd::sunset_propagator{

  template<int n_omega,
	   std::size_t l_x, std::size_t l_y, std::size_t n_t>
  
  fermi_G_t<l_x, l_y, n_t>
  
  project_to_frequency(fermi_G_t<l_x, l_y, n_t> const& Sigma){
    Real const Beta = 2.*Sigma[0].bounds_average;
    Real const two_over_Beta = 2./Beta;
    Real const omega_n = 2.*ffd::core_math::Pi*(.5+n_omega)/Beta;
    fermi_G_t<l_x, l_y, n_t> ret;

    
    for(std::size_t r=0; r<size(Sigma); ++r){
      [[maybe_unused]] auto [x, y] = ffd::math_tools::DivMod(r, (1ul<<l_x));
      auto const Sigma_r = Sigma[r];
      Real const ReSigma = two_over_Beta*ffd::gauss_kronrod::
      	Integrate_63([=](Real const t){return Sigma_r(t)*std::cos(omega_n*t);},
      		     {0., Beta});
      Real const ImSigma = two_over_Beta*ffd::gauss_kronrod::
	Integrate_63([=](Real const t){return Sigma_r(t)*std::sin(omega_n*t);},
		     {0., Beta});
      ret[r] = 
	fermi_cheby_s<n_t>([=](Real t){
			     return (ReSigma*std::cos(omega_n*t)+
				     ImSigma*std::sin(omega_n*t));
			   }, ffd::make_array(0., Beta));

    }


    return ret;
  }


}//namespace
