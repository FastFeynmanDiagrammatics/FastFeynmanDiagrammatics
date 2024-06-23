

namespace ffd::inverse_matrix{

  template<typename Field>
  void
  Row1Plus2Row3Matrix4(int row1_,
		       Field scalar_,
		       int row2_,
		       std::vector<Field>& M_){


    
    int const size = ffd::core_math::sqrt_int(std::size(M_));

    

    for(int j = 0; j < size; ++j){
      M_[j+row1_*size] += scalar_*M_[j+row2_*size];
    }

    
    
  }
  


}//namespace
