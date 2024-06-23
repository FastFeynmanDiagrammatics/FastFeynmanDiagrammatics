namespace ffd::user_space{
  
  template<int order,
	   phys::alpha_shift_type alpha_t,
	   typename alpha_arg_t,
	   typename matrix_t>
  void
  alpha_shift(matrix_t& M,
	      alpha_arg_t alpha){
    if constexpr(alpha_t == phys::magic){
	
      }

  }


}//namespace
