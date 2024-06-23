namespace ffd::chebyshev_fft{

  struct ChebyshevFFT{
    std::array<Real, 2> LowerUpperBound;
    std::vector<Real> Coef;

    

    ChebyshevFFT(Real scalar): Coef(1, scalar*Real(2)) { for(int j:{0,1}) LowerUpperBound[j] = -1+2*j; }
    ChebyshevFFT(): ChebyshevFFT(Real(0)) {}

    
    template<typename function_t>
    ChebyshevFFT(function_t const& FunctionToBeApproximated,
		 std::array<Real, 2> lower_upper_bound,
		 Real AbsolutePrecision = 1e-10,
		 int Log2StartingOrder = 2);

    
    
    
    ChebyshevFFT& operator=(ChebyshevFFT const&) = default;    

    

    
    Real operator()(Real x) const;


    Real ErrorEstimate() const{return std::max(std::abs(Coef[size(Coef)-1]),
					       std::abs(Coef[size(Coef)-1-(size(Coef)>1)]));}

		 

    ChebyshevFFT Derivative() const;

    ChebyshevFFT Derivative(int derivative_order) const;

    

    ChebyshevFFT Integral() const;
    
  };
  
  std::size_t size(ChebyshevFFT const& P){ return size(P.Coef);}
  
}//namespace
