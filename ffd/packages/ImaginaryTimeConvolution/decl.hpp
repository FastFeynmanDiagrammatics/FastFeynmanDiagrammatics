namespace ffd::imaginary_time_convolution{

  
  template<typename Field,
	   typename poly_t = ffd::chebyshev_polynomial::ChebyshevPolynomial<Field>>
  poly_t
  Convolute1And2FromZeroTo3(std::function<Field(Real)> f,
			    std::function<Field(Real)> g,
			    Real Beta_,
			    bool is_fermion_,
			    Real absolute_precision_ = 1e-10);

  
}//namespace

