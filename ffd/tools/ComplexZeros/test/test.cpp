#include"../../../std.hpp"
#include"../../../core/Math.hpp"
#include"../../AnalyticalDerivative.hpp"
#include"../ComplexZeros.hpp"


using namespace ffd::User;


template<typename F>
F Z(Real Beta, F mu, F U){
  return 1. + exp(mu*Beta)*2. + exp(mu*Beta*2.-Beta*U);
}


template<typename F>
F Density(Real Beta, F mu, F U){
  return (exp(mu*Beta)+exp(mu*Beta*2.-Beta*U))*2./Z(Beta, mu, U);
}



int main(){

  Real beta = 1;
  Real mu0 = -0.31;
  Real n0 = Density(beta, mu0, 0.); 
  
  auto func = [=](Complex U){return Z(beta, mu0+U*n0*.5, U);};

  ffd::ComplexZeros::ComplexZeros Func(func);

  std::cout<<Func.FindNearestZero()<<std::endl;
  

}



