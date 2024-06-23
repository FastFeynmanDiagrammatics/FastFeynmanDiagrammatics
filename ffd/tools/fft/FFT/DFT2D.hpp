

namespace ffd::fft{

  std::vector<Complex>

  DFT2D(std::vector<Complex> f,
	unsigned long size_x = 0ul){
    using ulong_t = unsigned long;

    
    ulong_t const size_f = size(f);
    if( size_x == 0ul ){
      size_x = ffd::core_math::sqrt_int( size_f );
    }
    std::vector<Complex> tilde_f(size_f, 0.);
    assert(( size_f%size_x == 0 ));
    ulong_t const size_y = size_f/size_x;
    

    auto phase_x = create_phases(size_x);
    auto phase_y = create_phases(size_y);
    

    for( ulong_t kx=0; kx < size_x; ++kx ){
      for( ulong_t ky=0; ky < size_y; ++ky ){
	for( ulong_t x=0; x < size_x; ++x ){
	  for( ulong_t y=0; y < size_y; ++y ){
	    auto x_kx = (x*kx)%size_x;
	    auto y_ky = (y*ky)%size_y;
	    tilde_f[kx + size_x*ky] +=
	      phase_x[x_kx]*phase_y[y_ky]*f[x+size_x*y];
	  }
	}
      }
    }

    
    return tilde_f;
  }

}//namespace
