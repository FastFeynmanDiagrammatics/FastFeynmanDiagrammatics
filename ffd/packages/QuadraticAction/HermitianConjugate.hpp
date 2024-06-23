

namespace ffd::user_space{

  template<typename Field>
  auto HermitianConjugate(QuadraticActionTerm<Field> const& A_){
    QuadraticActionTerm<Field> ret;
    ret = Bar(A_);
    ret.Value = ffd::complex_conjugate::ComplexConjugate(ret.Value);
    return ret;
  }

  
  template<typename Field>
  auto HermitianConjugate(QuadraticAction<Field> const& A_){
    QuadraticAction<Field> ret;
    for(std::size_t j=0; j < size(A_); ++j){
      ret.push_back(HermitianConjugate(A_[j]));
    }
    return ret;
  }


  
}//namespace
