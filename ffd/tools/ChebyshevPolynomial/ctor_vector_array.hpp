namespace ffd::chebyshev_polynomial{

  template<typename Field>
  ChebyshevPolynomial<Field>::ChebyshevPolynomial(std::vector<Field> values_at_nodes,
						  std::array<Real, 2> low_upp):
    ChebyshevPolynomial(values_at_nodes, low_upp[0], low_upp[1]) {}


}//namespace
