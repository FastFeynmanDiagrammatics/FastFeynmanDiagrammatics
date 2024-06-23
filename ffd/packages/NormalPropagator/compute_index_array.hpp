namespace ffd::normal_propagator{

  template<int d,
	   typename vector_t1,
	   typename vector_t2
	   >
  std::array<int, 2>
  compute_index_array(vector_t1 x,
		      vector_t2 L){
    int index = 0, product_L = 1;
    for(int j=0; j<d; ++j){
      index += product_L*x[j];
      product_L *= L[j];
    }
    return {index, product_L};
  }


}//namespace
