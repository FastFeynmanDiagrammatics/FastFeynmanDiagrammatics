

namespace ffd::inverse_matrix::unit_test{

  template<typename Field>
  std::vector<Field>
  invert_one_minus_A(std::vector<Field> A_){
    std::vector<Field> inverse(std::size(A_), 0.);


    
    std::vector<Field> A_pow_n = A_;


    
    int const size = ffd::core_math::sqrt_int(std::size(A_));


    
    for(int j=0; j < size; ++j){
      inverse[j+j*size] = 1.;
    }

    

    const int NumIterations = 1000;
    for(int k=0; k < NumIterations; ++k){
      for(int j=0; j < size*size; ++j){
	inverse[j] += A_pow_n[j];
      }
      
      A_pow_n = multiply_matrices(A_pow_n, A_);
    }


    
    return inverse;
  }


}//namespace
