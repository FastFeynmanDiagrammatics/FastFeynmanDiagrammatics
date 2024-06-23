namespace ffd::heat_bath_mc::unit_test{

  template<int order,
	   int order_plus,
	   int half_num_sites>

  
  void invoke_translation_invariant(Real ext_f = 1.,
				    unsigned long N_iter = 1ul<<18,
				    Real omega_f = 1.33){
    auto g = [=](Real x){
	       return std::cos(3.1415926535*omega_f*std::abs(x))/(1+ext_f*x);
	     };

    
      auto f =
	[=](auto r){
	  Real norm2 = 0.;
	  for(int j=0; j<size(r); ++j){
	    for(int k=0; k<size(r); ++k){
	      norm2 += std::pow(r[j]-r[k],2);
	    }
	  }
	  return g(norm2);
	    //g(arg_exp);
	    //std::exp(-std::sqrt(arg_exp))*sign;
	};
      
      exact_all_orders<order, order_plus, half_num_sites, true>(f, N_iter);
  }


}//namespace
