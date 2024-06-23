

namespace ffd::chebyshev_polynomial_2d::unit_test{
  
  void creation(){
    using namespace std;
    
    Real x_m =0, x_M = 1;
    Real y_m =0, y_M = 1;
    Real precision = 1e-10;
    int order = 40;

    auto f = [](Real x, Real y){return (exp(-x)+2*exp(-y))*sin(10*x*y+.3*sin(5*x +y));};
    
    ChebyshevPolynomial2D<Real> P(f,
				  order,
				  {{{x_m, x_M}, {y_m, y_M}}});

    ChebyshevPolynomial2D<Real> P2(f,
				   {{{x_m, x_M}, {y_m, y_M}}},
				   precision);

  }


}//namespace
