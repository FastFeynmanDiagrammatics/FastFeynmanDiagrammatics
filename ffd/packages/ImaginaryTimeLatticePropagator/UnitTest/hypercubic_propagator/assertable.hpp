

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  template<int d>
  bool assertable_hypercubic_propagator(std::array<int, d> L,
					std::array<Real, 2> Mu,
					std::array<std::array<Real, d>, 2> t_spin_j,
					Real beta){
    bool IsOk = true;


    HypercubicPropagatorExact<d> G0;
    NambuPropagatorExact<d> N0;

    G0.L = L;
    N0.L = L;
    G0.Mu = Mu;
    N0.Mu = Mu;
    G0.t_spin_j = t_spin_j;
    N0.t_spin_j = t_spin_j;
    G0.beta = beta;
    N0.beta = beta;
    N0.Delta = 0;


    ffd::lattice::HypercubicLattice<d> Hypercubic(L);
    ffd::imaginary_time::PeriodicImaginaryTime BetaStrip(beta);

    

    auto O = CreateCoordinates(BetaStrip,
			       Hypercubic);
    Real tau0=1e-15;
    component<0>(O) = tau0;

    
    QuadraticAction<Real> S0;
    for( int spin: {-1, 1} ){
      S0 += Bar(Psi_(spin))(O)*Psi_(spin)(O) * ( - Mu[(1+spin)/2] );
    }


    auto hopping = hypercubic_hopping<d, 0>(O, t_spin_j);
    S0 += hopping;

    
    
    auto G0_f = ImaginaryTimeLatticePropagator(BetaStrip,
					       Hypercubic,
					       S0);


    int N_samples_tau = 5;
    for(int r=0; r<4; ++r){
      for(int t=0; t<N_samples_tau; ++t){
	Real tau = (t+.2)*beta/N_samples_tau;
	std::array<int, d> rr;
	rr.fill(0);
	rr[0] = r;

	// auto psi = Nambu(Psi_(1)[0]);


	
	for(int spin: {-1, 1} ){
	  auto psi = Psi_(spin)[0];
	  auto X = O;
	  component<0>(X) += tau;
	  component<1>(X) += r;
	  auto pair_psi = std::pair{psi, std::any(X)};
	  auto pair_bpsi = std::pair{Bar(psi), std::any(O)};

	  // std::cerr<<spin<<" "<<r<<" "<<tau/beta<<" "<<
	  //   G0(tau, rr, spin)<<" "<<G0(true, tau, rr, spin)-G0(tau, rr, spin);
	  // std::cerr<<" "<<G0(true, tau, rr, spin)-spin*N0(spin*tau, rr, {{(1-spin)/2, (1-spin)/2}});
	  // std::cerr<<" "<<G0(true, tau, rr, spin)+G0_f({{pair_psi, pair_bpsi}})<<std::endl;
	  
	  IsOk = IsOk &&
	    std::abs(G0(true, tau, rr, spin) + G0_f({{pair_psi, pair_bpsi}})) <
	    20*std::numeric_limits<Real>::epsilon();
	  
	}
      }
    }
    

    return IsOk;
  }

};
