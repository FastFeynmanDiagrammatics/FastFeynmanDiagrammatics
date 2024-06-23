

namespace ffd::inverse_matrix::unit_test{

  template<typename Field>
  std::vector<Field>
  multiply_matrices(std::vector<Field> const& M1,
		    std::vector<Field> const& M2){
    std::vector<Field> product(std::size(M1), 0.);


    int const size = ffd::core_math::sqrt_int(std::size(M1));
    

    for(int j = 0; j < size; ++j){
      for(int k = 0; k < size; ++k){
	for(int m = 0; m < size; ++m){
	  product[j*size+k] += M1[j*size+m]*M2[m*size+k];
	}
      }
    }
    return product;
  }


}//namespace
