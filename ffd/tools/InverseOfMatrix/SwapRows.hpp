

namespace ffd::inverse_matrix{

  template<typename Field>
  void
  SwapRows(std::vector<Field>& M_,
	   std::array<int, 2> rows){

    int const size = ffd::core_math::sqrt_int(std::size(M_));

    
    for(int j = 0; j < size; ++j){
      std::swap(M_[j + rows[0]*size], M_[j + rows[1]*size]);
    }
    
  }


}//namespace
