namespace ffd::chebyshev_polynomial{

  template<typename Field>
  ChebyshevPolynomial<Field>::ChebyshevPolynomial(std::function<Field(Real)> function_,
						  int N_Cheby, 
						  std::array<Real, 2> lower_upper_limit):
    ChebyshevPolynomial<Field>(function_, N_Cheby, lower_upper_limit[0], lower_upper_limit[1]){}

}//namespace
