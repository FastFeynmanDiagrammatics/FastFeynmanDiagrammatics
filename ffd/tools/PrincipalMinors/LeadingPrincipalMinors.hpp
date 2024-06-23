namespace ffd::principal_minors{

  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  LeadingPrincipalMinors(matrix_t&& matrix){

    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
    using namespace ffd::nilpotent_polynomial;

    const BinaryInt order = ffd::core_math::sqrt_int(size(matrix));
    assert (order*order == (int)size(matrix));

    std::vector<element_t> minors(order+1, 0.);
  	minors[0] = 1.;
  	minors[1] = matrix[0];

  	for (BinaryInt l=0; l<order; ++l){
  		const element_t inv = 1./minors[l+1];
  		for (BinaryInt i = l+1; i<order; ++i){
  			const element_t inv1 = matrix[l+i*order] * inv;
  			for (BinaryInt j = l+1; j<order; ++j){
  				matrix[j+i*order] -= matrix[j+l*order] * inv1;
  			}
  		}
  		minors[l+2] = matrix[l+1+(l+1)*order];
  	}

    for (BinaryInt l = 1; l < order; ++l){
      minors[l+1] *= minors[l];
    }

    // std::cout << "n:" << std::endl;
    // for (BinaryInt level = 0; level <= order; ++level){
    //   std::cout << level << " " << minors[level] << std::endl;
    // }

    return minors;
  }

  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  LeadingPrincipalMinors2n(matrix_t&& matrix){

    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
    using namespace ffd::nilpotent_polynomial;

    const BinaryInt order = ffd::core_math::sqrt_int(size(matrix));
    assert (order*order == (int)size(matrix));

    // auto matrix = matrix_;

    std::vector<element_t> minors(order+1, 0.);
    minors[0] = 1.;
    minors[2] = matrix[0]*matrix[order+1] - matrix[1]*matrix[order];

    for (BinaryInt l=0; l<order; l+=2){
      auto inv_0 =  matrix[l+1+(l+1)*order];
      auto inv_1 = -matrix[l+1+l*order];
      auto inv_2 = -matrix[l+(l+1)*order];
      auto inv_3 =  matrix[l+l*order];
      const element_t det = inv_0*inv_3 - inv_1*inv_2;
      inv_0 /= det;
      inv_1 /= det;
      inv_2 /= det;
      inv_3 /= det;
      for (BinaryInt i = l+2; i<order; ++i){
        const element_t m1 = matrix[l+i*order];
        const element_t m2 = matrix[l+1+i*order];
        const element_t inv1_0 = m1*inv_0 + m2*inv_2;
        const element_t inv1_1 = m1*inv_1 + m2*inv_3;
        for (BinaryInt j = l+2; j<order; ++j){
          matrix[j+i*order] -= matrix[j+l*order]*inv1_0 + matrix[j+(l+1)*order]*inv1_1;
        }
      }
      minors[l+2] = det;
    }

    for (BinaryInt l = 2; l < order; ++l){
      minors[l+2] *= minors[l];
    }

    // std::cout << "2n:" << std::endl;
    // for (BinaryInt level = 0; level <= order; ++level){
    //   std::cout << level << " " << minors[level] << std::endl;
    // }
    // exit(1);

    return minors;
  }

}//namespace
