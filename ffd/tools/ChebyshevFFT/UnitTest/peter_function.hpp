

namespace ffd::chebyshev_fft::unit_test{

  struct peter_function{
    
    Real operator()(Real x) const{
      using std::sin, std::exp, std::pow;
      return exp(-pow(.1, 100)/pow(x, 100))*(2+sin( 100*std::pow(x, .2342342) + 97*pow(x, 10)))*
	sin(81*std::tanh(sin(100.*x)) + 73*sin(121.*sin(39*sin(73*x*x*x + 57*sqrt(x)))));
    }
    
  };


}//namespace
