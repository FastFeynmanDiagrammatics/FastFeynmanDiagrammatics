namespace ffd::chebyshev_fft{


  Real NormMax(ChebyshevFFT const& P_){
    Real norm = 0;
    for(std::size_t j=0; j<size(P_); ++j){
      norm = std::max(std::abs(P_.Coef[j]), norm);
    }
    return norm;
  }


}//namespace
