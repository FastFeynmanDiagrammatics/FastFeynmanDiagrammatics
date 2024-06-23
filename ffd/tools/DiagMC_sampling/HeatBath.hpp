namespace ffd::diagmc_sampling{

  template<order_t Order,
	   order_t OrderMinus,
	   typename polynomial_t,
	   typename seed_function_vector_t,
	   bool sum_order_n_minus_1>
	   
  std::tuple<vertex_number_t,
	     signed_time_low_order_t,
	     signed_time_high_order_t>
  
  HeatBath(polynomial_t P,
	   seed_function_vector_t seeds){
    using ffd::vector_range::Range;

    
    constexpr static unsigned long two_Order =
      1ul<<((unsigned)Order);
    

    if constexpr(OrderMinus == 1){
	std::array<Real, Order> stack_cumul;
	for( auto j: Range(Order) ){
	  stack_cumul[j]  =  ( (j==0) ? 0. : stack_cumul[j-1] );
	  stack_cumul[j]  +=  std::abs( P[two_Order-1-(1<<j)] )*seeds[j];
	}
	
	Real signed_time_low_order = 0.;
	Real stack_sum = stack_cumul[Order-1];
	for( auto j: Range(Order) ){
	  stack_cumul[j] /= stack_sum;
	  signed_time_low_order += std::copysign(stack_cumul[j], P[two_Order-1-(1<<j)]);
	}
	
	
	Real const signed_time_high_order = P[two_Order-1]/stack_sum;
	int const chosen_stack = SelectProbabilityStack(stack_cumul);
	if constexpr(!sum_order_n_minus_1){
	    signed_time_low_order = std::copysign(Real(1), P[ two_Order - 1 - (1<<chosen_stack) ]);
	}

	
	return std::make_tuple(chosen_stack,
			       signed_time_low_order,
			       signed_time_high_order);


	
      }else if constexpr(OrderMinus == 2){

	
	std::vector<Real> stack_cumul(Order*Order, 0.);
	for( auto j: Range(Order) ){
	  for( auto k: Range(Order) ){
	    int const jk = k+j*Order;
	    stack_cumul[jk]  =  ( (jk==0) ? 0. : stack_cumul[jk-1] );
	    if( j != k ){
	      std::size_t const subset = two_Order - 1 - (1<<j) - (1<<k);
	      stack_cumul[jk]  +=  std::abs( P[subset] )*seeds[jk];
	    }
	  }
	}

    
	Real const stack_sum = stack_cumul[ Order*Order - 1 ];
	for( auto jk: Range(Order*Order) ){
	  stack_cumul[jk] /= stack_sum;
	}
	

	
	Real const signed_time_high_order = P[two_Order-1]/stack_sum;
	int const chosen_stack = SelectProbabilityStack(stack_cumul);
	int const j_chosen = chosen_stack % Order;
	int const k_chosen = chosen_stack / Order;
	Real const P_low = P[two_Order-1-(1<<j_chosen)-(1<<k_chosen)];
	Real const signed_time_low_order = std::copysign(1., P_low);
	
	
	return std::make_tuple(chosen_stack,
			       signed_time_low_order,
			       signed_time_high_order);

	
      }else{
      assert(false);
      return std::make_tuple(0, 0., 0.);
    }
  }
  
}//namespace
