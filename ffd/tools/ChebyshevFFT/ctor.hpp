namespace ffd::chebyshev_fft{

  template<typename function_t>
  ChebyshevFFT::
  ChebyshevFFT(function_t const& F,
	       std::array<Real, 2> lower_upper_bound,
	       Real AbsolutePrecision,
	       int Log2StartingOrder):
    LowerUpperBound(lower_upper_bound){
    int log2_order = Log2StartingOrder;
    do{
      Coef = ChebyshevInterpolateFFT(F, LowerUpperBound, log2_order);
      ++log2_order;
    }while( ErrorEstimate() > AbsolutePrecision );
  }

}//namespace
