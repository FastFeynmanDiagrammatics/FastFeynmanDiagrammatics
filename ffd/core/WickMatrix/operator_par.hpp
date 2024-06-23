namespace ffd::wick_matrix{

  template<typename Field, typename intType>  
  Field WickMatrix<Field, intType>::operator()(intType j_, intType k_) const{
    if(j_ == k_ && !NotHermitian && IsFermion){
      return Field();
    }
    if(!NotHermitian && j_ > k_ ){
      if(IsFermion){
// 	if constexpr (std::is_same_v<Field, ffd::feynman_edge::FeynmanEdge> ){
// return (*this)(k_, j_);
// }else{
	return -(*this)(k_, j_);
      }else{
	return (*this)(k_, j_);
      }
    }
    return Components[ConvertToInternalRepresentation(j_, k_)];
  }

  
  template<typename Field, typename intType>
  Field& WickMatrix<Field, intType>::operator()(intType j_, intType k_, [[maybe_unused]] const char assign_){
    assert(assign_ == 'a');
    assert(NotHermitian || !IsFermion || j_ < k_);
    if(!NotHermitian && j_ > k_){
      return (*this)(k_, j_, 'a');
    }
    return Components[ConvertToInternalRepresentation(j_, k_)];
  }

  
}//namespace ffd::wick_matrix
