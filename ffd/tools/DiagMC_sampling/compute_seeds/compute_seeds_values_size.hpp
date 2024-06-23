namespace ffd::diagmc_sampling{


  template<typename seed_function_t,
	   typename origin_t,
	   typename vector_t1,
	   typename vector_t2,
	   typename int_t>

  void
  compute_seed_values_size(seed_function_t seed_,
			   origin_t O_,
			   vector_t1 X_,
			   vector_t2& ret,
			   int_t size_ret){
    using ffd::vector_range::Range;
    for(int_t j=0; j<size_ret; ){
      ret[j*size_ret] = 0;
      for( auto m: Range(size_ret) ){
	auto& ret_element = ret[j*size_ret+m];
	if( j != m ){
	  ret_element = seed_(X_[j], X_[m]);
	}else{
	  ret_element = seed_(X_[j], O_ );
	}
      }
    }
  }

}//namespace
