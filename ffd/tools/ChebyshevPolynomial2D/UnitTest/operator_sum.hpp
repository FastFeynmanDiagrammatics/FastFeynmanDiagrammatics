

namespace ffd::chebyshev_polynomial_2d::unit_test{

  void operator_sum(){
    using namespace std;

    
    Real a= 2.12312312, b=1.21312, c=.20923231;
    Real precision;
    


    auto f = [a, b](Real x, Real y){return sin(a*x*x + b*y*log(1+y));};
    auto g = [a, b, c](Real x, Real y){return cos(exp(a*x*x + b*y*log(c+y)));};


    
    auto F = ChebyshevPolynomial2D<Real>(f,
					 {{{0, 1}, {0, 2}}},
					 precision=1e-7);
    // std::cerr<<F.Order<<std::endl;

    
    
    auto G = ChebyshevPolynomial2D<Real>(g,
					 {{{0, 1}, {0, 2}}},
					 precision=1e-7);
    //    std::cerr<<G.Order<<std::endl;

    
    
    auto F_plus_G = F + G;
    // std::cerr<<F_plus_G.Order<<std::endl;


    
    F += G;
    

    
    int N_samples = 12;
    for(int y=0; y<N_samples; ++y){
      for(int x=0; x<N_samples; ++x){
	Real xx = (.1+x)/N_samples;
	Real yy = 2*(.1+y)/N_samples;
	//	std::cerr<<(F_plus_G(xx, yy) - f(xx, yy) - g(xx, yy))<<std::endl;
	assert(std::abs( F_plus_G(xx, yy) - f(xx, yy) - g(xx, yy) ) < precision);
	assert(std::abs( F(xx, yy) - f(xx, yy) - g(xx, yy) ) < precision);
      }
    }

  }

}//namespace
