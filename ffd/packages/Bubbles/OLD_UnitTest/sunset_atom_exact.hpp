

namespace ffd::rpa_ladder::unit_test{


  struct G_sunset_atom{

    Real beta, mu0, U;
    

    Real operator()(Real tau){
      Real G_sunset = 0.;

      
      Real S = std::abs(U)/(2*std::cosh(beta*mu0/2));
      for(int j: {-1, 1}){
	Real E = -mu0 + j*S;
	G_sunset += -.5*exp(-tau*E)/(1+exp(-beta*E));
      }
      return G_sunset;
    }

    
  };


}//namespace
