namespace ffd::chebyshev_polynomial_2d{

  template<typename Field>
  struct ChebyshevPolynomial2D{
    std::vector<Field> Coef;
    int Order;
    std::array<std::array<Real, 2>, 2> Domain;


    
    ChebyshevPolynomial2D(): Coef(1, 0.), Order(1) {}

    ChebyshevPolynomial2D(std::function<Field(Real, Real)> FunctionToBeApproximated,
			  int OrderOfThePolynomialApproximation, 
			  std::array<std::array<Real, 2>, 2> rectangular_domain);

    ChebyshevPolynomial2D(std::function<Field(Real, Real)> FunctionToBeApproximated,
			  std::array<std::array<Real, 2>, 2> rectangular_domain,
			  Real AbsoluteRequestedPrecision = 1e-10,
			  int MinimalAcceptedOrder = 2,
			  Real AddendIncreasingIterativeOrder = 1.,
			  Real FactorIncreasingIterativeOrder = 1.3);

    
    

    Field operator()(Real, Real) const;

    Real ErrorEstimate() const;
    
  };
	
	
}//namespace
