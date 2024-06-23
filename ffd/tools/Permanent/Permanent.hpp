namespace ffd::permanent{

  template<typename matrix_t>
   
  typename std::decay<decltype(std::declval<matrix_t>()[0])>::type
  Permanent(matrix_t&& matrix){

  using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
    
    using namespace ffd::nilpotent_polynomial;
    using ffd::vector_range::Range;

  	const int linear_size = ffd::core_math::sqrt_int(size(matrix));
    assert (linear_size*linear_size == (int)size(matrix));
    const int cardinality_2 = (1 << (linear_size-1));


    if (linear_size==1){
      return matrix[0];
    }

    element_t total = 0.0;
    std::vector<element_t> row_sum(linear_size);
    int old_gray = 0;
    element_t sign = 1.0;

    double row_sum_prod = sign;
    for (uint i = 0; i<linear_size; ++i){
      row_sum[i] = 0.0;
      for (uint j = 0; j<linear_size; ++j){
        row_sum[i] += matrix[i+linear_size*j];
      }
      row_sum_prod *= row_sum[i];
    }

    for (uint set=1; set<cardinality_2+1; ++set){
      total += row_sum_prod;
      int new_gray = set ^ (set/2);
      int gray_diff = old_gray ^ new_gray;

      int gray_diff_index = 0;
      if (gray_diff)
      {
        while (gray_diff>>=1) gray_diff_index+=1;
      } 

      double direction = (old_gray>new_gray ? 2.0 : -2.0);
      

      for (uint i = 0; i<linear_size; ++i){
        row_sum[i] += matrix[i+linear_size*gray_diff_index] * direction;
      }

      sign = -sign;
      old_gray = new_gray;
      
       //this loop could be acceletared for NilP's by fast summation
      row_sum_prod = sign;
    	for (uint i = 0; i<linear_size; ++i){
    		row_sum_prod *= row_sum[i];
    	}
    }
    return total/(double)cardinality_2;
  }
} //namespace
