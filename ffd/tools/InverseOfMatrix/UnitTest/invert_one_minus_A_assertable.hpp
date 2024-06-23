

namespace ffd::inverse_matrix::unit_test{

  template<typename Field>
  bool
  invert_one_minus_A_assertable(std::vector<Field> A_,
				Real precision = -.1){
    bool IsOk = true;


    
    int const size = ffd::core_math::sqrt_int(std::size(A_));

    
    
    Real Precision = precision > 0 ?
      precision :
      10*std::numeric_limits<Real>::epsilon();



    auto inverse1 = invert_one_minus_A(A_);


    
    auto one_minus_A = A_;
    for(int j=0; j < size; ++j){
      one_minus_A[j + j*size]  -=  1.;
    }
    
    
    
    for(int j=0; j < size*size; ++j){
      one_minus_A[j] *= -1.;
    }


    
    auto inverse2 = InverseOfMatrix(one_minus_A);

    

    for(int j = 0; j < size*size; ++j){
      IsOk = IsOk &&
      	std::abs(inverse1[j] - inverse2[j]) <
      	Precision;
      // std::cerr<<inverse1[j]<<" "<<inverse2[j]<<std::endl;
    }

    
    
    return IsOk;
  }

  
}//namespace
