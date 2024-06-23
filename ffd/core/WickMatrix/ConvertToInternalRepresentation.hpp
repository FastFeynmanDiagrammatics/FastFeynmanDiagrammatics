

namespace ffd::wick_matrix{

  template<typename Field, typename intType>
  intType WickMatrix<Field, intType>::ConvertToInternalRepresentation(intType j_, intType k_) const{
    intType const& n = this->NumberDestructionOperators;
    if(this->NotHermitian){
      return (j_ * n)  +  (k_);
    }else if(this->IsFermion){
      return (j_ * n - (j_*(j_+1))/2)  +  (k_ - (j_+1));
    }else{
      return (j_ * n - (j_*(j_-1))/2)  +  (k_ - j_);;
    }
  }


}//namespace ffd::wick_matrix
