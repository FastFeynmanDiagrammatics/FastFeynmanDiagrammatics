namespace ffd::user_space{
  
  //TODO	
  template<int order,
	   phys::alpha_shift_type alpha_t,
	   typename alpha_arg_t = phys::default_tmp_t,
	   typename matrix_t = phys::default_tmp_t>
  void
  alpha_shift(matrix_t& M,
	      alpha_arg_t alpha = phys::default_tmp_t());
  
}//namespace
