namespace ffd::principal_minors{

  template<uint max_order = 0, typename matrix_t>
  // template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  InversePrincipalMinors(matrix_t&& matrix, BinaryInt step=1){

    using namespace ffd::determinant;
    using namespace ffd::nilpotent_polynomial;
    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

    const int linear_size = ffd::core_math::sqrt_int(size(matrix)/(step*step));
    assert (linear_size*linear_size*step*step == (int)size(matrix));


    uint max_ord;
    if (max_order == 0){
      max_ord = linear_size;
    }
    else{
      max_ord = max_order;
    }

    const int cardinality = (1 << linear_size);

    std::vector<element_t> minors(cardinality);
    minors[0] = 1.0;

    for (uint set = 1; set < cardinality; ++set){
      // std::cout << set << std::endl;
      uint pop_set = __builtin_popcount(set);
      if (pop_set<=max_ord)
      {
        uint set_ = ffd::set_theory::BinarySpread(set, step);
        uint w_in = 0;
        uint minor_size = step*pop_set;
        std::vector<element_t> minor_matrix(minor_size*minor_size);
        for (uint v_in = 0; v_in < linear_size*step; ++v_in){
          int w_ou = 0;
          if (set_ & (1 << v_in)){
            for (uint v_ou = 0; v_ou < linear_size*step; ++v_ou){
              if (set_ & (1 << v_ou)){
                minor_matrix[w_in*minor_size + w_ou] = matrix[v_in*linear_size*step + v_ou];
                w_ou++;
              }
            }
            w_in++;
          }
        }
	auto [I, det] = ffd::inverse_matrix_gauss::
	  Inverse_and_Determinant(minor_matrix);
        minors[set] = det;
      }
    }

    return minors;
  }


}//namespace
