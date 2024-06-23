namespace ffd::imaginary_time_convolution_2d{

  //returns the correct answer
  //  O N L Y    I F 
  // tau1 <= tau2
  //otherwise returns the analytic continuation of tau1 <= tau2
  template<typename Field,
	   typename chebyshev_2d_t>

  chebyshev_2d_t

  FastFermionicConvolute123FromZeroTo4(std::function<Field(Real)> f_,
				       std::function<Field(Real)> g_,
				       std::function<Field(Real)> h_,
				       Real beta_,
				       Real absolute_precision_){
    auto first_integral =
      [f_, g_, h_, absolute_precision_](Real tau1, Real tau2){
	auto f_G_H_fixed12 = [tau1, tau2, f_, g_, h_](Real tau){
			       return f_(tau)*g_(tau1-tau)*h_(tau2-tau);
			     };
	return ffd::integrate_clenshaw_curtis::
	  IntegrateWithClenshawCurtis<Field>(f_G_H_fixed12,
					     {{0, tau1}},
					     absolute_precision_);
      };									   




    auto second_integral =
      [f_, g_, h_, absolute_precision_, beta_](Real tau1, Real tau2){
	auto f_G_H_fixed12 = [tau1, tau2, f_, g_, h_, beta_](Real tau){
			       return -f_(tau)*g_(beta_+tau1-tau)*h_(tau2-tau);
			     };
	return ffd::integrate_clenshaw_curtis::
	  IntegrateWithClenshawCurtis<Field>(f_G_H_fixed12,
					     {{tau1, tau2}},
					     absolute_precision_);
      };									   



    auto third_integral =
      [f_, g_, h_, absolute_precision_, beta_](Real tau1, Real tau2){
	auto f_G_H_fixed12 = [tau1, tau2, f_, g_, h_, beta_](Real tau){
			       return f_(tau)*g_(beta_+tau1-tau)*h_(beta_+tau2-tau);
			     };
	return ffd::integrate_clenshaw_curtis::
	  IntegrateWithClenshawCurtis<Field>(f_G_H_fixed12,
					     {{tau2, beta_}},
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
