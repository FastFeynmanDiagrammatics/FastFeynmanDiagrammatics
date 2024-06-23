

namespace ffd::chebyshev_fft::unit_test{

  void
  hard_fit(){

    Real precision = 1e-12;
    hard_function f;

    
    ChebyshevFFT P(f,
		   {-1., 2.},
		   precision);

    
    // std::cerr<<size(P)<<std::endl;
    int N_samples = 1000;
    for( int j=0; j<N_samples; ++j){
      Real x = -1 + (j+.13242342312312)*2./N_samples;
      assert((
	      std::abs( f(x)-P(x) ) < 10*precision
	      ));
      // std::cerr<<x<<" "<<std::setprecision(15)<<P(x)<<" "<<f(x)-P(x)<<std::endl;
    }
		   

  }

}//namespace
