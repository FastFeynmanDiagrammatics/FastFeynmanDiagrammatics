

namespace ffd::chebyshev_fft{

  ChebyshevFFT ChebyshevFFT::Derivative() const{
    ChebyshevFFT derivative = *this;
    auto const size_this = size(*this);
    
    using ffd::vector_range::Range;
    for( auto j: Range(size_this) ){
      derivative.Coef[j] = 0;
    }


    derivative.Coef[size_this-1] = 0.;
    derivative.Coef[size_this-2] = 2.*(size_this - 1)*Coef[size_this-1];
    for( auto j = size_this - 2; j > 0; --j ){
      derivative.Coef[j-1] = derivative.Coef[j+1] + 2.*j*Coef[j];
    }
    Real diff_interval = return_half_mean_diff(LowerUpperBound)[1];
    for( auto j: Range(size_this) ){
      derivative.Coef[j] /= diff_interval;
    }

    
    return derivative;
  }


}//namespace
