

namespace ffd::chebyshev_2d_fft::unit_test{

  void
  hard_fit(){
    Real precision = 1e-12;

    hard_function_x f_x;

    Chebyshev2DFFT P(f_x,
		     std::array{std::array{0., 1.},
				  std::array{0., 1.}},
		     precision);

    
    auto [size_x, size_y] = sizes(P);
    // std::cerr<<"sizes_p = "<<size_x<<" "<<size_y<<std::endl;

    int const N_samples = 20;
    for(int x=0; x<N_samples; ++x){
      for(int y=0; y<N_samples; ++y){
	Real xx = (x+.12312312)/N_samples;
	Real yy = (y+.12312312)/N_samples;
	assert((
		std::abs( P(xx, yy) - f_x(xx, yy) ) <
		10*precision
		));
	// std::cerr<<P(xx, yy)<<" "<<P(xx, yy) - f_x(xx, yy)<<std::endl;
      }
    }
    
  }


}//namespace
