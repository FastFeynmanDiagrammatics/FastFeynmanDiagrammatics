namespace ffd::chebyshev_2d_fft{
  
  struct Chebyshev2DFFT{
    using ulong_t = unsigned long;
    inline static const Real AbsolutePrecisionDefault = 1e-10;
    inline static const ulong_t Log2StartingSizeDefault = 2;
    ulong_t Log2XSize;
    std::array<std::array<Real, 2>, 2> LowerUpperBounds;
    std::vector<Real> Coef;
    ulong_t Log2YSize() const{ assert((size(Coef)%(1ul<<Log2XSize)==0));
	return ffd::math_tools::log2_int(size(Coef)/(1ul<<Log2XSize));}


    
    Chebyshev2DFFT(Real scalar): Log2XSize(0ul), Coef(1, scalar){
      for(int j:{0,1}) for(int u:{0,1}) LowerUpperBounds[j][u] = -1+2*u;}
    Chebyshev2DFFT(): Chebyshev2DFFT(0.) {}

    template<typename function2d_t>
    Chebyshev2DFFT(function2d_t Function2DToBeApproximated,
		   decltype(LowerUpperBounds) lower_upper_bounds,
		   Real AbsolutePrecision = AbsolutePrecisionDefault,
		   ulong_t Log2StartingSize = Log2StartingSizeDefault);

    // template<typename function2d_t>
    // Chebyshev2DFFT(function2d_t Function2DToBeApproximated,
    // 		   std::array<Real, 2> lower_upper_bound,
    // 		   Real AbsolutePrecision = AbsolutePrecisionDefault,
    // 		   ulong_t Log2StartingSize = Log2StartingSizeDefault):
    //   Chebyshev2DFFT(Function2DToBeApproximated,
    // 		     std::array{lower_upper_bound, lower_upper_bound},
    // 		     AbsolutePrecision,
    // 		     Log2StartingSize) {}

    

    Real operator()(std::array<Real, 2> X) const;
    Real operator()(Real x, Real y) const{ return this->operator()({x, y}); }
    
    
    
    std::array<Real, 2>
    ErrorEstimateXY();
    
    Real ErrorEstimate(){
      auto const& [e_x, e_y] = ErrorEstimateXY();
      return std::max({e_x, e_y, std::abs( Coef[ size(Coef) - 1])});
    }

  };


  std::size_t size(Chebyshev2DFFT const& P){
    return size(P.Coef);
  }

  std::array<std::size_t,2>
  sizes(Chebyshev2DFFT const& P){
    return std::array{1ul<<P.Log2XSize,
			1ul<<P.Log2YSize()};
  }
  
	
}//namespace
