

namespace ffd::chebyshev_2d_fft{

  std::array<Real, 2>

  Chebyshev2DFFT::
  ErrorEstimateXY(){
    std::array<Real, 2> xy_errors;
    xy_errors.fill(0.);

    
    using std::abs;
    

    ulong_t const log2ysize = Log2YSize();
    ulong_t const pow_2_n = 1<<Log2XSize;
    ulong_t const pow_2_m = 1<<log2ysize;

    
    for( ulong_t kx = 0; kx < pow_2_n-1; ++kx){
      xy_errors[1] = std::max({xy_errors[1],
			       abs(Coef[kx+pow_2_n*(pow_2_m-2)]),
			       abs(Coef[kx+pow_2_n*(pow_2_m-1)])});
    }

    
    for( ulong_t ky = 0; ky < pow_2_m-1; ++ky){
      xy_errors[0] = std::max({xy_errors[0],
			       abs(Coef[pow_2_n-2+pow_2_n*ky]),
			       abs(Coef[pow_2_n-1+pow_2_n*ky])});
    }

    
    return xy_errors;
  }


}//namespace
