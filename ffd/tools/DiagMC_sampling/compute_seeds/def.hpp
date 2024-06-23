namespace ffd::diagmc_sampling{
  
  template<typename seed_function_t,
	   typename origin_t,
	   typename vector_t>
    
  auto
  compute_seed_values(seed_function_t seed_,
		      origin_t O_,
		      vector_t X_){
    std::size_t const size_X = size( X_ );
    std::vector<Real> ret(size_X*size_X);
    

    compute_seed_values_size(seed_, O_, X_, ret, size_X);


    return ret;
  }


  
  template<int Order,
	   typename seed_function_t,
	   typename origin_t,
	   typename vector_t>
    
  auto
  compute_seed_values(seed_function_t seed_,
		      origin_t O_,
		      vector_t X_){
    std::array<Real, Order*Order> ret;

 
    compute_seed_values_size(seed_, O_, X_, ret, Order);
    
    
    return ret;
  }

  
}//namespace
