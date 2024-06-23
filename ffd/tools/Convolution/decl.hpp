namespace ffd::convolution{


  namespace parameters{
    Real const AbsolutePrecision = 1e-10;
  }

  
  template<typename function1_t,
	   typename function2_t>

  ffd::chebyshev_fft::ChebyshevFFT
  
  Convolution(function1_t f, function2_t g,
	      std::array<Real, 2> interval,
	      Real AbsolutePrecision = parameters::AbsolutePrecision,
	      bool IsBosonic = true){
    using chebyshev_t =  ffd::chebyshev_fft::ChebyshevFFT;


    Real interval_range = interval[1]- interval[0];
    

    chebyshev_t fc(f, interval, AbsolutePrecision);
    chebyshev_t gc(g, interval, AbsolutePrecision);


    auto fg = [fc, gc](Real x, Real y)->Real{ return fc( x - y )*gc( y );};
    
    
    auto conv_fg_first_integral =
      [interval, fg, AbsolutePrecision](Real x)->Real{
	auto fg_fixedx = [x, fg](Real y)->Real{return fg(x, y);};
	return ffd::integrate_clenshaw_curtis::
	  IntegrateWithClenshawCurtis<Real>(fg_fixedx,
					    {interval[0], x-interval[0]},
					    AbsolutePrecision);
      };
    

    
    Real zeta_sign = 2*(IsBosonic) - 1;
    auto conv_fg_second_integral =
      [=](Real x)->Real{
	auto fg_fixedx = [=](Real y)
			 {return zeta_sign*fg(interval_range+x-y, y);};
	return ffd::integrate_clenshaw_curtis::
	  IntegrateWithClenshawCurtis<Real>(fg_fixedx,
					    {x-interval[0], interval[1]},
					    AbsolutePrecision);
      };


    

    auto conv_fg = [conv_fg_first_integral, conv_fg_second_integral](Real x)->Real{
		      return conv_fg_first_integral(x)+
			conv_fg_second_integral(x);
		    };


    
    return ffd::chebyshev_fft::ChebyshevFFT(conv_fg, interval, AbsolutePrecision);

  }
	
}//namespace
