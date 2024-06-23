namespace ffd::integrate_clenshaw_curtis{

  template<typename Field = Real>
  Field IntegrateWithClenshawCurtis(std::function<Field(Real)> FunctionToIntegrate,
				    std::array<Real, 2> IntervalOfIntegration,
				    Real RequestedAbsolutePrecision = 1e-10,
				    int MinimalOrder = 4,
				    Real AddendIncreasingIterativeOrder = 1.,
				    Real FactorIncreasingIterativeOrder = 1.51){
    int iterative_order = MinimalOrder;
    Real iterative_order_Real = iterative_order;
    Real interval = IntervalOfIntegration[1]-IntervalOfIntegration[0];
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Field>
      poly(FunctionToIntegrate, iterative_order, IntervalOfIntegration);
    while(RequestedAbsolutePrecision < 2*poly.ErrorEstimate()*std::abs(interval)/pow(size(poly), 2)){
      poly = ffd::chebyshev_polynomial::ChebyshevPolynomial<Field>
	(FunctionToIntegrate, iterative_order, IntervalOfIntegration);
      iterative_order_Real *= FactorIncreasingIterativeOrder;
      iterative_order_Real += AddendIncreasingIterativeOrder;
      iterative_order = std::round(iterative_order_Real);
    }
    Real ret = poly.Coef[0];
    for(int j=1; 2*j < size(poly); ++j){
      ret -= 2*poly.Coef[2*j]/((2*j+1)*(2*j-1));
    }
    return .5*ret*interval;
  }

}//namespace
