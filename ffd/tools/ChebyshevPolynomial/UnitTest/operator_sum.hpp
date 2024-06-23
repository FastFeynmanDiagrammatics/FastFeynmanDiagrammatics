

namespace ffd::chebyshev_polynomial::unit_test{

  void operator_sum(){
    using namespace std;
    
    Real a = -2.131231231;
    Real b = 5.598234;
    Real precision = 1e-10;
    
    auto f = [](Real x){return exp(-sin(x))*cos(x);};
    auto g = [](Real x){return cos(sin(x))/(1+x*x);};
    auto f_plus_g = [f, g](Real x){return f(x)+g(x);};

    ChebyshevPolynomial<Real> F(f, {a, b}, precision);
    ChebyshevPolynomial<Real> G(g, {a, b}, precision);
    ChebyshevPolynomial<Real> F_PLUS_G(f_plus_g, {a, b}, precision);

    
    auto F_plus_G = F + G;

    
    const int NumIterations = 100;
    for(int j=0; j < NumIterations; ++j){
      Real tau = .5*(a+b) + .5*j*(b-a)/NumIterations;
      assert( abs( F_PLUS_G(tau) - F_plus_G(tau) ) < precision );
    }

    
  }

}//namespace
