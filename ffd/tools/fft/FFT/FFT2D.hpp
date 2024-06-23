namespace ffd::fft{

  std::vector<Complex>

  FFT2D(std::vector<Complex> f,
	unsigned long size_x){
    using ulong_t = unsigned long;

    
    ulong_t const size_f = size(f);
    if( size_x == 0ul ){
      size_x = ffd::core_math::sqrt_int( size_f );
    }
    std::vector<Complex> f_kx_ky(size_f), f_kx_y(size_f);
    assert(( size_f%size_x == 0 ));
    ulong_t const size_y = size_f/size_x;

    

    for( ulong_t y = 0; y < size_y; ++y ){
      std::vector<Complex> f_x(size_x);
      int const shift_y = y*size_x;
      for( ulong_t x=0; x < size_x; ++x ){
	f_x[x] = f[x + shift_y];
      }

      
      auto fourier_x = FFT(f_x);
      for( ulong_t kx=0; kx < size_x; ++kx ){
	f_kx_y[kx+shift_y] = fourier_x[kx];
      }
    }


    
    for( ulong_t kx=0; kx < size_x; ++kx ){
      std::vector<Complex> f_y(size_y);
      for( ulong_t y=0; y < size_y; ++y ){
	f_y[y] = f_kx_y[kx + y*size_x];
      }
      
      
      auto fourier_y = FFT(f_y);
      for( ulong_t ky=0; ky < size_y; ++ky ){
	f_kx_ky[kx + size_x*ky] = fourier_y[ky];
      }
    }
    

    return f_kx_ky;
  }


}//namespace
