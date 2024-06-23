namespace ffd::imaginary_time_convolution{

  template<typename Field,
	   typename poly_t>
  poly_t
  Convolute1And2FromZeroTo3(std::function<Field(Real)> f_,
			    std::function<Field(Real)> g_,
			    Real beta_,
			    bool is_fermion_,
			    Real absolute_precision_){

    auto f_g = [f_, g_](Real tau_1, Real tau_2)->Real{ return f_(tau_1 - tau_2)*g_(tau_2);};


    
    auto conv_f_g_first_integral =
      [f_g, absolute_precision_](Real tau_1)->Real{
	auto f_g_fixed1 = [tau_1, f_g](Real tau_2)->Real{return f_g(tau_1, tau_2);};
	return ffd::integrate_clenshaw_curtis::IntegrateWithClenshawCurtis<Real>(f_g_fixed1,
										 {0, tau_1},
										 absolute_precision_);
      };


    
    Real zeta_sign = 1-2*(is_fermion_);
    auto conv_f_g_second_integral =
      [beta_, zeta_sign, f_g, absolute_precision_](Real tau_1)->Real{
	auto f_g_fixed1 = [tau_1, zeta_sign, beta_, f_g](Real tau_2)->Real{return zeta_sign*f_g(beta_+tau_1, tau_2);};
	return ffd::integrate_clenshaw_curtis::IntegrateWithClenshawCurtis<Real>(f_g_fixed1,
										 {tau_1, beta_},
										 absolute_precision_);
      };


    

    auto conv_f_g = [conv_f_g_first_integral, conv_f_g_second_integral](Real tau_1)->Real{
		      return conv_f_g_first_integral(tau_1)+
			conv_f_g_second_integral(tau_1);
		    };


    
    poly_t ret(conv_f_g, {0, beta_}, absolute_precision_);
    return ret;
  }


}//namespace
