

namespace ffd::rpa_ladder::unit_test{
  using namespace std;

  
  struct AtomicLobster{

    Real Beta, mu0, U;

    
    AtomicLobster(){}

    
    inline Real ReturnA() const{
      Real z = std::exp(Beta*mu0);
      return (1-z*z)/std::pow(1+z, 2);
    }


    
    inline Real ReturnD(Real tau1, Real tau2) const{
      Real A = ReturnA();
      return -U*exp((tau1+tau2)*mu0)/(pow(1+exp(Beta*mu0), 2)*(1-exp(-Beta*(A*U- 2*mu0))));
    }



    Real operator()(Real tau1, Real tau2) const{
      Real D = ReturnD(tau1, tau2);
      Real A = ReturnA();
      Real second_term =
	( 1 + exp(Beta*mu0) )*( exp(Beta*mu0-tau2*A*U) - exp(-tau1*A*U) );
      Real third_term = -exp(Beta*(2*mu0-A*U));
      return D*( 1 + second_term + third_term );
    }
    
    
  };

}//namespace
