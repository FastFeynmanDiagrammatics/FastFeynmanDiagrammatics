namespace ffd::chebyshev_fft{

  ChebyshevFFT
  operator+(ChebyshevFFT const& P1_,
	    ChebyshevFFT const& P2_){
    using std::size;
    
    ChebyshevFFT const& P_small_ = size(P1_) < size(P2_) ? P1_ : P2_;
    ChebyshevFFT const& P_big_ = size(P1_) < size(P2_) ? P2_ : P1_;
    ChebyshevFFT sum = P_big_;
    for(std::size_t j=0; j < size(P_small_); ++j){
      sum.Coef[j] += P_small_.Coef[j];
    }
    return sum;
  }


  

  auto
  operator-(ChebyshevFFT const& P1_,
	    ChebyshevFFT const& P2_){
    return P1_ + P2_*(-1.);
  }


  
  ChebyshevFFT& operator+=(ChebyshevFFT      & P1_,
			   ChebyshevFFT const& P2_){
    P1_ = P1_ + P2_;
    return P1_;
  }



}//namespace
