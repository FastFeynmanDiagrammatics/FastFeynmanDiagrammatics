namespace ffd::chebyshev_polynomial{

  template<typename Field>
  auto
  operator*(ChebyshevPolynomial<Field> const& P_,
	    Field scalar_){
    ChebyshevPolynomial<Field> product = P_;
    
    for(int j=0; j < size(P_); ++j){
      product.Coef[j] *= scalar_;
    }
    return product;
  }


}//namespace
