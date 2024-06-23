namespace ffd::chebyshev_polynomial{

  template<typename Field>
  ChebyshevPolynomial<Field>::ChebyshevPolynomial(std::function<Field(Real)> function_,
						  int N_Cheby, 
						  Real lower_limit,
						  Real upper_limit){
    auto x_Cheby = ReturnsChebyshevNodesOfOrder1From2To3(N_Cheby, lower_limit, upper_limit);
    std::vector<Field> vector_of_values(N_Cheby);
    for(int j=0; j < N_Cheby; ++j){
      vector_of_values[j] = function_(x_Cheby[j]);
    }
    *this = ChebyshevPolynomial(vector_of_values, lower_limit, upper_limit);
  }


}//namespace
