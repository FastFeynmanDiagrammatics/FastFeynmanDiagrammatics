namespace ffd::imaginary_time_convolution_2d{

  //returns the correct answer
  //  O N L Y    I F 
  // tau1 <= tau2
  //otherwise returns the analytic continuation of tau1 <= tau2
  template<typename chebyshev_2d_t,
	   typename f_t,
	   typename g_t,
	   typename h_t>	   

  chebyshev_2d_t

  VeryFastFermionicConvolute123FromZeroTo4(f_t&& f_v,
					   g_t&& g_v,
					   h_t&& h_v,
					   Real beta,
					   Real precision){
    auto first_integral =
      [&](Real tau1, Real tau2){
	auto integrand = [&](Real tau){
	  return f_v(tau)*g_v(tau1-tau)*h_v(tau2-tau);
	};
	return ffd::gauss_kronrod::
	Integrate(integrand, {0, tau1}, precision);
      };									   



    auto second_integral =
      [&](Real tau1, Real tau2){
	auto integrand = [&](Real tau){
	  return -f_v(tau)*g_v(beta+tau1-tau)*h_v(tau2-tau);
	};
	return ffd::gauss_kronrod::
	Integrate(integrand, {tau1, tau2}, precision);
      };									   


    
    auto third_integral =
      [&](Real tau1, Real tau2){
      auto integrand = [&](Real tau){
	return f_v(tau)*g_v(beta+tau1-tau)*h_v(beta+tau2-tau);
      };
      return ffd::gauss_kronrod::
      Integrate(integrand, {tau2, beta}, precision);
    };									   


    
    auto integral_function = [&](Real tau1, Real tau2){
      return first_integral(tau1, tau2)+
      second_integral(tau1, tau2)+
      third_integral(tau1, tau2);
    };

    

    return chebyshev_2d_t(integral_function,
			  {{{0., beta}, {0., beta}}},
			  precision);

  }


}//namespace
