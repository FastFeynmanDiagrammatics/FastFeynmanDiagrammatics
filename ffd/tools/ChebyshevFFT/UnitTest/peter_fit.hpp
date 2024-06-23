

namespace ffd::chebyshev_fft::unit_test{

  void
  peter_fit(){
    peter_function f;
    Real precision = 1e-12;

    
    ChebyshevFFT P(f,
		   {0., 1.},
		   precision);

    
    // ffd::chebyshev_polynomial::
    //   ChebyshevPolynomial<Real> P2(f,
    // 				   {0., 1.},
    // 				   precision);
    

    // std::cerr<<size(P)<<std::endl;
    int N_samples = 1000;
    for( int j=0; j<N_samples; ++j){
      Real x = (j+.13242342312312)*1./N_samples;
      assert((
	      std::abs( f(x)-P(x) ) < 10*precision
	      ));
      // std::cerr<<x<<" "<<std::setprecision(15)<<P(x)<<" "<<f(x)-P(x)<<std::endl;
    }
		   

    
  }


}//namespace
