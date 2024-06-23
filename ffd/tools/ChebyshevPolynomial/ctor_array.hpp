namespace ffd::chebyshev_polynomial{

  template<typename Field>
  template<std::size_t order>
  
  ChebyshevPolynomial<Field>::ChebyshevPolynomial(std::array<Field, order> const& values_at_nodes,
						  std::array<Real, 2> low_upp):
    Coef(size(values_at_nodes), 0.),
    LowerLimit(low_upp[0]),
    UpperLimit(low_upp[1]){


    using ffd::core_math::Pi;


    std::size_t const size_Coef = size(Coef);
    
    for(std::size_t j=0; j < size_Coef; j++){
      for(std::size_t k=0; k < size_Coef; k++){
	Coef[j] += values_at_nodes[k]*cos(Pi*j*(k+0.5)/size_Coef);
      }
      Coef[j] *= 2./size_Coef;
    }
  }

}//namespace
