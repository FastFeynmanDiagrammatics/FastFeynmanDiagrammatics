namespace ffd::bubbles::unit_test{

  bool gamma0_atom(Real beta, Real mu0, Real U, [[maybe_unused]] Real alpha = .7, Real precision = 1e-9){
    bool IsOk = true;
    
    using chebyshev_t = ffd::chebyshev_fft::ChebyshevFFT;
    
    using namespace std;


    int NumIterations = 200;
    auto g0 = [beta, mu0](Real tau)->Real{return -exp(tau*mu0)/(1+exp(beta*mu0));};
    chebyshev_t g0_p(g0, {0, beta}, precision);
    std::array<chebyshev_t, 1> G0;
    G0[0] = g0_p;
    std::vector<decltype(G0)> G0_vec;
    G0_vec.push_back(G0);
    std::array<decltype(G0_vec), 2> G0_spin;
    G0_spin.fill(G0_vec);
    
    
    
    ffd::imaginary_time::PeriodicImaginaryTime IG(beta);
    ffd::lattice::HypercubicLattice<1> OneSite({1});
    

    auto x = ffd::lattice::GetLinearSizes(OneSite);
    
    Bubbles
      Bub(IG,
	  OneSite,
	  G0_spin,
	  U);

    
    for(int j=0; j < NumIterations; ++j){
      [[maybe_unused]] auto error_iteration = Bub.IterateP();
      // std::cerr<<"Max_error(iteration="<<j<<") = "<<error_iteration<<std::endl;
      error_iteration = 0;
    }

    
    
    Real A_bubbles = (1-exp(2*beta*mu0))/pow(1+exp(beta*mu0), 2);
    auto P = [A_bubbles, U, mu0, beta](Real tau)->Real{
		    return -A_bubbles*exp(-tau*(A_bubbles*U-2*mu0))/(1-exp(-beta*(A_bubbles*U-2*mu0)));};



    
    const int N_grid = 100;
    for(int j=0; j<N_grid; ++j){
      Real tau = (j+.5)*beta/N_grid;
      // std::cerr<<"P("<<tau<<") = "<<P(tau)<<", Error2 =  "<<Bub.P[0][0](tau)-P(tau)<<std::endl;
      IsOk = IsOk && abs( Bub.P[0][0](tau)-P(tau) ) < precision;
      // std::cerr<<"P("<<tau<<") = "<<P(tau)<<", Error1 =  "<<Bubbles1.P[0](tau)-P(tau)<<std::endl;
      // IsOk = IsOk && abs( Bub.P[0](tau)-P(tau) ) < 2*precision;
    }
    
    return IsOk;
  }


}//namespace
 
