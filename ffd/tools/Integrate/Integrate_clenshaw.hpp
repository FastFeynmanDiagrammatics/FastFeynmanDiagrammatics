namespace ffd::integrate{

  template<typename function_t>
  Real
  Integrate_clenshaw(function_t Integrand,
		     std::array<Real, 2> interval,
		     Real RequestedAbsolutePrecision){

    
    Real const len_interval = interval[1] - interval[0];
    ffd::chebyshev_fft::ChebyshevFFT F(Integrand,
				       interval,
				       RequestedAbsolutePrecision*
				       len_interval);


    
    Real integral = F.Coef[0];
    for(std::size_t j=1; 2*j < size(F); ++j){
      integral -= 2.*F.Coef[2*j]/((2*j+1)*(2*j-1));
    }
    
    
    return .5*integral*len_interval;
  }
  

}//namespace
