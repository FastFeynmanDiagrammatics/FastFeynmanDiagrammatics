//ChebyshevPolynomial builds a polynomial, stored in coefficients (coef)
//You need to give the value of a function in particular points (see below)
//you call the Chebyshev approximation with operator()
//you can get an error estimate with ErrorEstimate()
namespace ffd::chebyshev_polynomial{

  template<typename Field = Real>
  class ChebyshevPolynomial{
  public:
    std::vector<Field> Coef;
    Real LowerLimit, UpperLimit;
    

    explicit ChebyshevPolynomial(Field x): Coef(1, x*Real(2)), LowerLimit(-1.), UpperLimit(1.) {}
    explicit ChebyshevPolynomial(): ChebyshevPolynomial(0.) {}

    //it builds a Chebyshev polynomial approximation 
    //of a function f(x), x\in[a, b] (default values are a=-1, b=1).
    //We need to give values at Chebyshev nodes
    //values_at_nodes[j] = f(.5*(a+b) + .5*(b-a)*cos(Pi*(j+0.5)/N_Cheby))
    //where j\in\{0,1,...,N_Cheby-1\}
    //In alternative, in order to determine the
    //nodes, we can use the free function
    //ReturnsChebyshevNodesOfOrder1From2To3(int, Real, Real)
    ChebyshevPolynomial(std::vector<Field> values_at_the_Chebyshev_nodes,
			Real lower_limit = -1.,
			Real upper_limit =  1.);
    ChebyshevPolynomial(std::vector<Field> values_at_the_Chebyshev_nodes,
			std::array<Real, 2> lower_upper_limit);
    template<std::size_t order>
    ChebyshevPolynomial(std::array<Field, order> const& values_at_the_Chebyshev_nodes,
			std::array<Real, 2> lower_upper_limit);
    
    
    //it builds a Chebyshev polynomial directly from the function,
    //order is given
    ChebyshevPolynomial(std::function<Field(Real)> FunctionToBeApproximated,
			int OrderOfThePolynomialApproximation, 
			Real lower_limit = -1.,
			Real upper_limit =  1.);
    ChebyshevPolynomial(std::function<Field(Real)> FunctionToBeApproximated,
			int OrderOfThePolynomialApproximation, 
			std::array<Real, 2> lower_upper_limit);

    
    //it builds a Chebyshev polynomial directly from the function,
    //order is determined by the sought error
    ChebyshevPolynomial(std::function<Field(Real)> FunctionToBeApproximated,
			Real lower_limit = -1.,
			Real upper_limit =  1.,
			Real AbsoluteRequestedPrecision = 1e-10,
			int MinimalAcceptedOrder = 4,
			Real AddendIncreasingIterativeOrder = 1.,
			Real FactorIncreasingIterativeOrder = 1.51);
    ChebyshevPolynomial(std::function<Field(Real)> FunctionToBeApproximated,
			std::array<Real, 2> lower_upper_limit,
			Real AbsoluteRequestedPrecision = 1e-10, 
			int MinimalAcceptedOrder = 4,
			Real AddendIncreasingIterativeOrder = 1.,
			Real FactorIncreasingIterativeOrder = 1.51);


    
      
    //it returns the Chebyshev-Polynomial approximation of your function
    //in a point x\in [a,b]
    Field operator()(Real x) const;

    
    //it returns an estimate of the maximal absolute error commited
    //in the interval [-1,1]
    Real ErrorEstimate() const{return std::max(std::abs(Coef[std::size(Coef)-1]),
					       std::abs(Coef[std::size(Coef)-1-(std::size(Coef)>1)]));}

    ChebyshevPolynomial& operator=(ChebyshevPolynomial const&) = default;

    
    template<typename T>
    friend std::size_t size(ChebyshevPolynomial<T> const&);

    
    void Purify(Real AbsolutePrecision = 1e-10);
    
    
  };

  
  template<typename T>
  std::size_t size(ChebyshevPolynomial<T> const& P){ return std::size(P.Coef);}

  
  
}//namespace ffd::chebyshev_polynomial


