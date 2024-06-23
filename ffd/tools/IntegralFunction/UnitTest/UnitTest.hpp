

namespace ffd::integral_function::unit_test{

  void UnitTest(){
    using namespace std;
    
    Real a = 2.12312312;
    Real gamma = 3.12312412;
    IntegralFunction<Real> F([a, gamma](Real x)->Real{return a*pow(x, gamma);}, 0);

    Real y = 0.543212312;
    assert( abs( F(y) - a*pow(y, gamma+1)/(gamma+1) ) < F.RequestedAbsolutePrecision );
    
    
    
    auto Sin = [a](Real x)->Real{return cos(a*x);};
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Real> P(Sin, {0, 1});
    IntegralFunction<Real> F0(P, 0);
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Real> Q(F0, {0, 1});

    assert( abs( Q(y) - sin(a*y)/a ) < F0.RequestedAbsolutePrecision );
    
    
    
    IntegralFunction<Real> F2([](Real x)->Real{return 2*x*sin(1/x)-cos(1/x);}, 0);
    F2.RequestedAbsolutePrecision = 1e-11;
    auto F2_ex = [](Real x)->Real{return x*x*sin(1/x);};

    assert( abs(F2(y)-F2_ex(y)) < 1e-3); //it is really struggling

  }
  
}//namespace
