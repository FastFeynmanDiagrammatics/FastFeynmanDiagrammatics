

namespace ffd::inverse_matrix{

  template<typename Field>
  void
  MultiplyRow1OfMatrix2By3(int row_,
			   std::vector<Field>& M_,
			   Field scalar_){


    
    int const size = ffd::core_math::sqrt_int(std::size(M_));



    for(int j = 0; j < size; ++j){
      M_[j+row_*size] *= scalar_;
    }
    

  }

}//namespace
