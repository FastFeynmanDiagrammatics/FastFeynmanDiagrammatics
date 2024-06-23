namespace ffd::diagmc_sampling{


  template<int Order, int OrderMinus,
	   typename vector_t>

  auto
  compute_pseudo_potentials(vector_t seeds_){
    if constexpr(OrderMinus == 1){
	std::array<Real, Order> ret;


	for( int j = 0; j < Order; ++j){
	  ret[j] = 0;
	  for( int m = 0; m < Order; ++m ){
	    ret[j] += seeds_[j*Order+m];
	  }
	  ret[j] /= Order;
	}

	
	return ret;
      }else if constexpr(OrderMinus == 2){
	std::array<Real, Order*Order> ret;


	auto pseudo1 = compute_pseudo_potentials<Order, 1>(seeds_);


	using ffd::vector_range::Range;
	for( int j: Range(Order) ){
	  for( int k: Range(Order) ){
	    if( j != k ){
	      Real pk = pseudo1[k];
	      Real pj = pseudo1[j];
	      ret[j*Order + k] = (Order/(Order-1.))*pj*pk-(.5/(Order-1.))*(pj+pk)*seeds_[Order*j+k];
	    }
	  }
	}
	

	return ret;
      }else{
      return;
    }
  }



}//namespace
