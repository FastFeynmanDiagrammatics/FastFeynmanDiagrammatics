namespace ffd::chebyshev_fft{

  ChebyshevFFT
  ChebyshevFFT::Integral() const{
    ChebyshevFFT integral = *this;
    auto const size_this = size(*this);
    using ffd::vector_range::Range;
    

    Real sign_jp1 = 1., integ = 0.;
    Real half_diff = .5*return_half_mean_diff(LowerUpperBound)[1];
    for(  auto  j:  Range( 1, size_this - 1 )  ){
      integral.Coef[j] = half_diff*( Coef[j-1] - Coef[j+1] )/j;
      integ += sign_jp1*integral.Coef[j];
      sign_jp1 = -sign_jp1;
    }
    integral.Coef[ size_this - 1 ]  =
      half_diff * Coef[ size_this - 2 ] / ( size_this - 1 );

    integ += sign_jp1 * integral.Coef[ size_this - 1 ];
    integral.Coef[0] = 2.*integ;
    
    
    return integral;
  }


}//namespace
