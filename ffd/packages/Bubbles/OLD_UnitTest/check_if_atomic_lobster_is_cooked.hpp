

namespace ffd::rpa_ladder::unit_test{

  template<bool IsLadder = true, bool EliminatingTadpoles = false>
  bool check_if_atomic_lobster_is_cooked(Real Beta,
					 Real mu0,
					 Real U,
					 Real alpha = .5,
					 Real precision = 1e-5){
    bool IsOk = true;

    auto G0 = [Beta, mu0](Real tau){return -exp(tau*mu0)/(1+exp(Beta*mu0));};
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Real> G0_element{G0, {{0, Beta}}, precision};
    std::vector<ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>> G0_vec;
    G0_vec.push_back(G0_element);
    


    RPA_Ladder<IsLadder, EliminatingTadpoles> Ladder(G0_vec,
						     Beta,
						     U,
						     alpha,
						     precision);

    

    Ladder.IterateUntilRequestedPrecision();



    
    auto Lobster = Ladder.ReturnLobster();



    AtomicLobster Lobster_exact;
    Lobster_exact.Beta = Beta;
    Lobster_exact.U = U;
    Lobster_exact.mu0 = mu0;

    

    int N_samples = 1000;
    for(int j1 = 0; j1 < N_samples; ++j1){
      for(int j2 = j1; j2 < N_samples; ++j2){
	Real t1 = (.2+j1)*Beta/N_samples;
	Real t2 = (.2+j2)*Beta/N_samples;
	// std::cerr<<Lobster_exact(t1, t2)<<" "<<Lobster_exact(t1, t2)-Lobster[0](t1, t2)<<std::endl;
	IsOk = IsOk && abs( Lobster_exact(t1, t2) - Lobster[0](t1, t2) )  <  precision;
      }
    }

    
    //    std::cerr<<Lobster[0].Order<<std::endl;

    
    return IsOk;
  }

}//namespace
