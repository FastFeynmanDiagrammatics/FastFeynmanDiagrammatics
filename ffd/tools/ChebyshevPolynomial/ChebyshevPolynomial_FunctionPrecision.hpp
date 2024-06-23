

namespace ffd::chebyshev_polynomial{

  template<typename Field>
  ChebyshevPolynomial<Field>::
  ChebyshevPolynomial(std::function<Field(Real)> func,
		      Real lower_limit,
		      Real upper_limit,
		      Real AbsolutePrecision, 
		      int MinimalAcceptedOrder,
		      Real AddendIncreasingIterativeOrder,
		      Real FactorIncreasingIterativeOrder){
    int iterative_order = MinimalAcceptedOrder;
    Real iterative_order_real = iterative_order;
    do{
      *this = ChebyshevPolynomial<Field>(func, iterative_order, lower_limit, upper_limit);
      iterative_order_real *= FactorIncreasingIterativeOrder;
      iterative_order = iterative_order_real + .5 + AddendIncreasingIterativeOrder;
    }while(AbsolutePrecision < ErrorEstimate() );
  }
    
      

}//namespace
