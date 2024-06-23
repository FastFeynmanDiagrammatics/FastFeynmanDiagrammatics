namespace ffd::fft{

  template<std::size_t nx,
	   std::size_t ny,
	   direction d = direct,
	   template<typename, std::size_t> typename array_t = std::array>
  
  array_t<Complex, (1ul<<(nx+ny))>

  FFT2D_l(array_t<Complex, (1ul<<(nx+ny))> const& f){
    static_assert( d == direct || d == inverse );
    std::size_t constexpr size_x_y = (1ul<<(nx+ny));
    std::size_t constexpr size_x = (1ul<<nx);
    std::size_t constexpr size_y = (1ul<<ny);

    
    array_t<Complex, size_x_y> f_kx_ky;
    array_t<Complex, size_x_y> f_kx_y;
    

    array_t<Complex, size_x> const phase_x = create_phases<size_x, d, array_t>();
    array_t<Complex, size_x> f_x, fourier_x;
    for(std::size_t y = 0; y<size_y; ++y){
      std::size_t const shift_y = y*size_x;
      for(std::size_t x=0; x<size_x; ++x){
	f_x[x] = f[x + shift_y];
      }

      
      fourier_x = FFT_l<nx, d, array_t>(f_x, phase_x);
      for(std::size_t kx=0; kx<size_x; ++kx){
	f_kx_y[kx+shift_y] = fourier_x[kx];
      }
    }


    array_t<Complex, size_y> const phase_y = create_phases<size_y, d, array_t>();
    array_t<Complex, size_y> f_y, fourier_y;
    for(std::size_t kx=0; kx<size_x; ++kx){
      for(std::size_t y=0; y<size_y; ++y){
	f_y[y] = f_kx_y[kx + y*size_x];
      }
      
      
      fourier_y = FFT_l<ny, d, array_t>(f_y, phase_y);
      for(std::size_t ky=0; ky<size_y; ++ky){
	f_kx_ky[kx + size_x*ky] = fourier_y[ky];
      }
    }
    
    
    return f_kx_ky;
  }



    template<std::size_t nx,
	     std::size_t ny,
	     direction d>
  
    std::array<Complex, (1ul<<(nx+ny))>

    FFT2D(std::array<Complex, (1ul<<(nx+ny))> const& f){
      return FFT2D_l<nx, ny, d, std::array>(f);
    }


}//namespace
