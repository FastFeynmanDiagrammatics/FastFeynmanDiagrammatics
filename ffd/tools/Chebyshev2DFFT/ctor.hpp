

namespace ffd::chebyshev_2d_fft{

  template<typename function2d_t>
  Chebyshev2DFFT::
  Chebyshev2DFFT(function2d_t F,
		 std::array<std::array<Real, 2>, 2> low_upp_bounds,
		 Real AbsolutePrecision,
		 unsigned long Log2StartingSize):
    LowerUpperBounds(low_upp_bounds){
    std::array<ulong_t, 2> log2_sizes;
    log2_sizes.fill(Log2StartingSize);
    do{
      Log2XSize = log2_sizes[0];
      Coef = ChebyshevInterpolate2DFFT(F, LowerUpperBounds, log2_sizes);
      auto [error_x, error_y] = ErrorEstimateXY();
      if( error_x < error_y ){
	++log2_sizes[1];
      }else{
	++log2_sizes[0];
      }
    }while( ErrorEstimate() > AbsolutePrecision );
  }

}//namespace
