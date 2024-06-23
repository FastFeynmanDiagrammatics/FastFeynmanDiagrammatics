

namespace ffd::rpa_ladder::unit_test{

  template<bool IsLadder = true, bool EliminatingTadpoles = false>
  bool sunset_atom(Real beta, Real mu0, Real U, Real alpha = .7, Real precision = 1e-9){
    bool IsOk = true;
    using namespace ffd::chebyshev_polynomial;

    
    auto G0 = [beta, mu0](Real tau){
		return -exp(tau*mu0)/(1+exp(beta*mu0));
	      };
    std::vector<ffd::chebyshev_polynomial::ChebyshevPolynomial<Real>> G0_vec;
    G0_vec.push_back( ChebyshevPolynomial<Real>(G0, {0, beta}, precision) );

    
    

    RPA_Ladder<IsLadder, EliminatingTadpoles> G_sunset(G0_vec,
						       beta,
						       U,
						       alpha,
						       precision);
    G_sunset.use_bold_G = false;


    auto Bubble = [G0](Real tau){return -G0(tau)*G0(tau);};


    G_sunset.P[0] = ChebyshevPolynomial<Real>(Bubble, {0, beta}, precision);
    
    
    const int N_iter = 200;
    for(int j=0; j < N_iter; ++j){
      G_sunset.IterateG();
    }

    
    
    G_sunset_atom G_exact{beta, mu0, U};


    const int N_samples = 30;
    for(int j=0; j < N_samples; ++j){
      Real tau = (j+.2)*beta/N_samples;
      //      std::cerr<<tau<<" "<<G0(tau)<<" "<<G_exact(tau)<<" "<<G_sunset.G[0](tau)<<std::endl;
      IsOk = IsOk && std::abs(G_exact(tau) - G_sunset.G[0](tau)) < precision;
    }
    return IsOk;
  }


}//namespace
