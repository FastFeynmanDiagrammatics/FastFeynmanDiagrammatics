namespace ffd::chebyshev_polynomial{

  template<typename Field>
  ChebyshevPolynomial<Field>::ChebyshevPolynomial(std::vector<Field> values_at_nodes,
						  Real lower_limit,
						  Real upper_limit): Coef(std::size(values_at_nodes), 0.),
								      LowerLimit(lower_limit),
								      UpperLimit(upper_limit){
    for(std::size_t j=0; j < std::size(Coef); j++){
      for(std::size_t k=0; k < std::size(Coef); k++){
	Coef[j] += values_at_nodes[k]*cos(ffd::core_math::Pi*j*(k+0.5)/std::size(Coef));
      }
      Coef[j] *= 2./std::size(Coef);
    }
  }


}//namespace
