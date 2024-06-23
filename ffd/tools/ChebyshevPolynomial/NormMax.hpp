namespace ffd::chebyshev_polynomial{

  template<typename Field>
  Real NormMax(ChebyshevPolynomial<Field> const& P_){
    Real norm = 0;
    for(int j=0; j<size(P_); ++j){
      norm = std::max(std::abs(P_.Coef[j]), norm);
    }
    return norm;
  }

}//namespace
