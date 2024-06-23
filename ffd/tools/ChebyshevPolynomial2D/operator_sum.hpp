namespace ffd::chebyshev_polynomial_2d{

  template<typename Field>
  ChebyshevPolynomial2D<Field>
  operator+(ChebyshevPolynomial2D<Field> const& addend1_,
	    ChebyshevPolynomial2D<Field> const& addend2_){
    bool first_addend_has_greater_order = addend1_.Order >= addend2_.Order;
    auto const& addend_max = first_addend_has_greater_order ? addend1_ : addend2_;
    auto const& addend_min = first_addend_has_greater_order ? addend2_ : addend1_;
    ChebyshevPolynomial2D<Field> sum = addend_max;
    for(int ky = 0; ky < addend_min.Order; ++ky){
      for(int kx = 0; kx < addend_min.Order; ++kx){
	sum.Coef[kx+ky*sum.Order] += addend_min.Coef[kx + ky*addend_min.Order];
      }
    }
    return sum;
  }



  template<typename Field>
  ChebyshevPolynomial2D<Field>&
  operator+=(ChebyshevPolynomial2D<Field>      & caller_,
	     ChebyshevPolynomial2D<Field> const& addend_){
    caller_ = caller_ + addend_;
    return caller_;
  }

  
}//namespace
