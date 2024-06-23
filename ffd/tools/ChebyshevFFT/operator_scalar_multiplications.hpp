namespace ffd::chebyshev_fft{

  auto
  operator*(ChebyshevFFT const& P_,
	    Real scalar_){
    ChebyshevFFT product = P_;
    
    for(std::size_t j=0; j < size(P_); ++j){
      product.Coef[j] *= scalar_;
    }
    return product;
  }


}//namespace
