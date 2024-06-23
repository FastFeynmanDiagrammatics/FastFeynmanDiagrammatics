

namespace ffd::rpa_ladder::unit_test{

  template<bool IsLadder = true, bool EliminatingTadpoles = false>
  bool check_if_lobster_is_cooked(Real beta,
				  Real mu0,
				  Real U,
				  int L,
				  Real alpha = .5,
				  Real precision = 1e-5){
    bool IsOk = true;

    
    int x_G0=0, y_G0=0;
    auto G0_tau = [mu0, L, beta, &x_G0, &y_G0](Real tau){return G0_square_lattice(x_G0,
										  y_G0,
										  tau,
										  beta,
										  mu0,
										  L);};

    
    
    std::vector<ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>> G0(L*L);
    for(y_G0=0; y_G0<L; ++y_G0){
      for(x_G0=0; x_G0<L; ++x_G0){
	G0[x_G0+y_G0*L] =
	  ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>(G0_tau,
							       {0, beta},
							       precision);
      }
    }
    
    
    RPA_Ladder<IsLadder, EliminatingTadpoles> Bubbles_auto(G0, beta, U, alpha, precision);
    Bubbles_auto.IterateUntilRequestedPrecision();

    //    std::cerr<<"start cooking lobster now..."<<std::endl;
    
    auto Lobster = Bubbles_auto.ReturnLobster();
    
    //    std::cerr<<Lobster[2](.2, .3)<<std::endl;
    
    return IsOk;
  }


}//namespace
