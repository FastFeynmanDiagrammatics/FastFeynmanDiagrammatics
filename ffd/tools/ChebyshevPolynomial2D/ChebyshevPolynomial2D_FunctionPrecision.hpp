

namespace ffd::chebyshev_polynomial_2d{

  template<typename Field>
  ChebyshevPolynomial2D<Field>::ChebyshevPolynomial2D(std::function<Field(Real, Real)>
						      FunctionToBeApproximated,
						      
  						      std::array<std::array<Real, 2>, 2>
						      rectangular_domain,
						      
  						      Real AbsoluteRequestedPrecision,
  						      int MinimalAcceptedOrder,
  						      Real AddendIncreasingIterativeOrder,
  						      Real FactorIncreasingIterativeOrder){
    int iterative_order = MinimalAcceptedOrder;
    Real iterative_order_real = iterative_order;
    do{
      *this = ChebyshevPolynomial2D(FunctionToBeApproximated,
				    iterative_order,
				    rectangular_domain);
      iterative_order_real *= FactorIncreasingIterativeOrder;
      iterative_order_real += AddendIncreasingIterativeOrder;
      iterative_order = .5 + iterative_order_real;
    }while( (*this).ErrorEstimate() > AbsoluteRequestedPrecision);
  }



}//namespace
