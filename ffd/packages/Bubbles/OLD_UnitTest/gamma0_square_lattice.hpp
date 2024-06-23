

namespace ffd::rpa_ladder::unit_test{

  template<bool IsLadder, bool EliminatingTadpoles>
  void gamma0_square_lattice(Real beta,
			     Real mu0,
			     Real U,
			     int L,
			     Real alpha,
			     Real precision){


    using namespace std;

    

    std::vector<ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>> G0(L*L);
    for(int y=0; y<L; ++y){
      for(int x=0; x<L; ++x){
	auto G0_tau =
	  [mu0, L, beta, x, y](Real tau){
	    return G0_square_lattice(x,
				     y,
				     tau,
				     beta,
				     mu0,
				     L);
	  };

	
	G0[x+y*L] =
	  ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>(G0_tau,
							       {0, beta},
							       precision);
      }
    }



    
    RPA_Ladder<IsLadder, EliminatingTadpoles> Bubbles_auto(G0, beta, U, alpha, precision);
    Bubbles_auto.IterateUntilRequestedPrecision();

    RPA_Ladder<IsLadder, EliminatingTadpoles> Bubbles(G0, beta, U, alpha, precision);

    const int NumIterations = 10;
    for(int j=0; j < NumIterations; ++j){
      auto error_iteration = Bubbles.IterateP(), error_iteration_G = 0.;
      if constexpr(EliminatingTadpoles){
	  error_iteration_G = Bubbles.IterateG();
	}
      // std::cerr<<"Max_error(iteration="<<j<<") = "<<error_iteration<<", "<<error_iteration_G<<std::endl;
      error_iteration = error_iteration_G;
    }
  }

}//namespace
