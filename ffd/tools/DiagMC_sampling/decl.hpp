namespace ffd::diagmc_sampling{

  using is_accepted_t = bool;
  using order_t = int;
  using vertex_number_t = int;
  using signed_time_low_order_t = Real;
  using signed_time_high_order_t = Real;

  
  
  template<order_t Order,
	   order_t OrderMinus,
	   typename polynomial_t,
	   typename seed_function_vector_t,
	   bool sum_order_n_minus_1 = true>
	   
  std::tuple<vertex_number_t,
	     signed_time_low_order_t,
	     signed_time_high_order_t>
  
  HeatBath(polynomial_t polynomial_filled_up_to_order_max,
	   seed_function_vector_t vector_of_seed_function_values);
	
	
}//namespace
