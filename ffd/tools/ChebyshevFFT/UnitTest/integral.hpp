

namespace ffd::chebyshev_fft::unit_test{

  void integral(){

    
    Real precision = 1e-12;

    Real a = -2., b = 3.;

    
    auto f = [](auto x){
	       using std::cos, std::sin;
	       return 2.*cos(x*x)-sin(x*x)/x/x;
	     };


    auto IIf = [](auto x){
	       using std::cos, std::sin;
	       return sin(x*x)/x;
	     };


    auto If = [IIf, a](auto x){
		return IIf(x) - IIf(a);
	      };
    

    ChebyshevFFT F(f,
		   {-2., 3.},
		   precision);
    
    
    auto IF = F.Integral();


    int N_samples = 400;
    using ffd::vector_range::Range;
    for( auto j: Range(N_samples) ){
      Real x = a + (j+.1)*(b-a)/N_samples;
      // std::cerr<<IF(x)<< " "<<IF(x)-If(x)<<std::endl;
      assert((
	      std::abs( IF(x) - If(x) ) <
	      precision
	      ));
    }
    
    
  }

  
}//namespace
