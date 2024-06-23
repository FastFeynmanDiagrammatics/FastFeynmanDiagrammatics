namespace ffd::heat_bath_mc::unit_test{

  template<int order,
	   int order_plus,
	   int half_num_sites>

  
  void invoke_exact_all_orders(Real ext_f = 1.,
			       unsigned long N_iter = 1ul<<18){
      auto f =
	[=](auto r){
	  Real arg_exp = 0.;
	  int sign = 1;
	  for(int j=0; j<size(r); ++j){
	    arg_exp += r[j]*r[j];
	    if( (r[j]&1) == 1){
	      sign = - sign;
	    }
	  }
	  return sign/(1.+ext_f*arg_exp);
	  //std::exp(-std::sqrt(arg_exp))*sign;
	};
      
      exact_all_orders<order, order_plus, half_num_sites, false>(f, N_iter);
  }

}//namespace
