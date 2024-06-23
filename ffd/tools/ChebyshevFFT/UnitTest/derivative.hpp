

namespace ffd::chebyshev_fft::unit_test{

  
  void derivative(){

    
    Real precision = 1e-12;

    
    auto f = [](Real x){
	       using std::exp;
	       return exp(-x*x);
	     };

    
    auto fp = [f](Real x){
		return -f(x)*x*2.;
	     };


    auto fpp = [f, fp](Real x){
		 return -f(x)*2. -2*x*fp(x);
	     };

    

    ChebyshevFFT P(f,
		   {-1., 2.},
		   precision);

    
    auto Pd = P.Derivative();

    auto Pdd = P.Derivative(2);
    

    ChebyshevFFT Pp(fp,
		    {-1., 2.},
		    precision);

    
    int const N_samples = 400;
    using ffd::vector_range::Range;
    for( auto j: Range(N_samples) ){
      auto x = -1. + (j+.2)*3./N_samples;
      assert((
	      std::abs(Pd(x) - fp(x)) < precision
	      ));
      assert((
	      std::abs(Pdd(x) - fpp(x)) < 10*precision
	      ));

      // std::cerr<<x<<" "<<Pd(x)<<" "<<Pd(x)-fp(x)<<std::endl;
      // std::cerr<<x<<" "<<Pdd(x)<<" "<<Pdd(x)-fpp(x)<<std::endl;
      
    }
    
  }


}//namespace
