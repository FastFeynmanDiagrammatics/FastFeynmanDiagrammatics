namespace ffd::chebyshev_polynomial{

  template<typename Field>
  ChebyshevPolynomial<Field>
  operator+(ChebyshevPolynomial<Field> const& P1_,
	    ChebyshevPolynomial<Field> const& P2_){
    using std::size;
    
    ChebyshevPolynomial<Field> const& P_small_ = size(P1_) < size(P2_) ? P1_ : P2_;
    ChebyshevPolynomial<Field> const& P_big_ = size(P1_) < size(P2_) ? P2_ : P1_;
    ChebyshevPolynomial<Field> sum = P_big_;
    for(int j=0; j < size(P_small_); ++j){
      sum.Coef[j] += P_small_.Coef[j];
    }
    return sum;
  }


  
  template<typename Field>
  auto
  operator-(ChebyshevPolynomial<Field> const& P1_,
	    ChebyshevPolynomial<Field> const& P2_){
    return P1_ + P2_*(-1.);
  }


  
  template<typename Field>
  ChebyshevPolynomial<Field>& operator+=(ChebyshevPolynomial<Field>      & P1_,
					 ChebyshevPolynomial<Field> const& P2_){
    P1_ = P1_ + P2_;
    return P1_;
  }
  
  

}//namespace
