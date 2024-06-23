

namespace ffd::rpa_ladder::unit_test{

  template<bool IsLadder = true, bool EliminatingTadpoles = false>
  bool gamma0_atom(Real beta, Real mu0, Real U, Real alpha = .7, Real precision = 1e-9){
    bool IsOk = true;

    
    using namespace std;


    int NumIterations = 200;
    auto g0 = [beta, mu0](Real tau)->Real{return -exp(tau*mu0)/(1+exp(beta*mu0));};
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Real> g0_p(g0, {0, beta}, precision);
    std::vector<ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>> G0;
    G0.push_back(g0_p);


    
    RPA_Ladder<IsLadder,
	       EliminatingTadpoles> Bubbles1(G0, beta, U, alpha, precision);
    Bubbles1.IterateUntilRequestedPrecision();
    
    
    
    RPA_Ladder<IsLadder,
	       EliminatingTadpoles> Bubbles2(G0, beta, U, alpha, precision);
    
    
    
    for(int j=0; j < NumIterations; ++j){
      auto error_iteration = Bubbles2.IterateP();
      //      std::cerr<<"Max_error(iteration="<<j<<") = "<<error_iteration<<std::endl;
      error_iteration = 0;
    }

    
    
    Real A_bubbles = (1-exp(2*beta*mu0))/pow(1+exp(beta*mu0), 2);
    auto P = [A_bubbles, U, mu0, beta](Real tau)->Real{
		    return -A_bubbles*exp(-tau*(A_bubbles*U-2*mu0))/(1-exp(-beta*(A_bubbles*U-2*mu0)));};



    
    const int N_grid = 100;
    for(int j=0; j<N_grid; ++j){
      Real tau = (j+.5)*beta/N_grid;
      //      std::cerr<<"P("<<tau<<") = "<<P(tau)<<", Error2 =  "<<Bubbles2.P[0](tau)-P(tau)<<std::endl;
      IsOk = IsOk && abs( Bubbles2.P[0](tau)-P(tau) ) < precision;
      //      std::cerr<<"P("<<tau<<") = "<<P(tau)<<", Error1 =  "<<Bubbles1.P[0](tau)-P(tau)<<std::endl;
      IsOk = IsOk && abs( Bubbles1.P[0](tau)-P(tau) ) < 2*precision;
    }
    
    return IsOk;
  }


}//namespace
