namespace ffd::imaginary_time_convolution_2d{

  template<typename Field,
	   typename chebyshev_2d_t>
  
  chebyshev_2d_t
  
  Convolute123FromZeroTo4(std::function<Field(Real)> f_,
			  std::function<Field(Real)> g_,
			  std::function<Field(Real)> h_,
			  Real beta_,
			  std::array<bool, 2> g_is_fermion__h_is_fermion_,
			  Real absolute_precision_){

    
    Real zeta_g = g_is_fermion__h_is_fermion_[0] == true ? -1. : 1.;
    Real zeta_h = g_is_fermion__h_is_fermion_[1] == true ? -1. : 1.;
    

    
    auto G_ = [g_, beta_, zeta_g](Real tau){
		return tau > 0 ? g_(tau) : zeta_g*g_(beta_+tau);
	      };
    auto H_ = [h_, beta_, zeta_h](Real tau){
		return tau > 0 ? h_(tau) : zeta_h*h_(beta_+tau);
	      };

    
    
    auto f_G_H_ = [f_, G_, H_](Real tau, Real tau1, Real tau2){
		    return f_(tau)*G_(tau1-tau)*H_(tau2-tau);
		  };


    

    auto first_integral =
      [f_G_H_, absolute_precision_](Real tau1, Real tau2){
	auto f_G_H_fixed12 = [tau1, tau2, f_G_H_](Real tau){
			       return f_G_H_(tau, tau1, tau2);
			     };
	return ffd::integrate_clenshaw_curtis::IntegrateWithClenshawCurtis<Field>(f_G_H_fixed12,
										  {{0, std::min(tau1, tau2)}},
										  absolute_precision_);
      };									   




    auto second_integral =
      [f_G_H_, absolute_precision_](Real tau1, Real tau2){
	auto f_G_H_fixed12 = [tau1, tau2, f_G_H_](Real tau){
			       return f_G_H_(tau, tau1, tau2);
			     };
	return ffd::integrate_clenshaw_curtis::IntegrateWithClenshawCurtis<Field>(f_G_H_fixed12,
										  {{std::min(tau1, tau2),
										    std::max(tau1, tau2)}},
										  absolute_precision_);
      };									   



    auto third_integral =
      [f_G_H_, absolute_precision_, beta_](Real tau1, Real tau2){
	auto f_G_H_fixed12 = [tau1, tau2, f_G_H_](Real tau){
			       return f_G_H_(tau, tau1, tau2);
			     };
	return ffd::integrate_clenshaw_curtis::IntegrateWithClenshawCurtis<Field>(f_G_H_fixed12,
										  {{std::max(tau1, tau2),
										    beta_}},
										  absolute_precision_);
      };									   




    auto integral_function = [first_integral, second_integral, third_integral](Real tau1, Real tau2){
			       return first_integral(tau1, tau2)+
				 second_integral(tau1, tau2)+
				 third_integral(tau1, tau2);
			     };




    return chebyshev_2d_t(integral_function,
			  {{{0, beta_}, {0, beta_}}},
			  absolute_precision_);
    
  }

}//namespace
