

namespace ffd::chebyshev_polynomial::unit_test{

  struct harder_than_not_so_easy{
    Real a = 10, b=10, c=3, d=1, e=3, f=-1;

    harder_than_not_so_easy() {}
    
    Real operator()(Real x) const{
      using std::sin, std::exp, std::abs;
      return sin(a*sin(b*sin(c*x+d)))*exp(e*x*x+f*x)/exp(e+abs(f));
    }
    
  };

}//namespace
