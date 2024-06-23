

namespace ffd::chebyshev_polynomial::unit_test{

  void operator_scalar_multiplication(){
    using namespace std;

    Real a = -.23123124213;
    Real b = 2.241867323;
    Real lambda = .123182132141;
    Real precision = 1e-10;

    auto f = [](Real x){return exp(-sin(x))*cos(x);};
    auto f_lambda = [f, lambda](Real x){return lambda*f(x);};

    ChebyshevPolynomial<Real> F(f, {a, b}, precision);
    ChebyshevPolynomial<Real> F_LAMBDA(f_lambda, {a, b}, precision);

    
    auto F_lambda = F*lambda;

    
    const int NumIterations = 100;
    for(int j=0; j < NumIterations; ++j){
      Real tau = .5*(a+b) + .5*j*(b-a)/NumIterations;
      assert( abs( F_LAMBDA(tau) - F_lambda(tau) ) < precision );
    }

    
  }

}//namespace
