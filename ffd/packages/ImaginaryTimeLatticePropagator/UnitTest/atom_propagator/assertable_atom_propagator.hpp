

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  auto
  assertable_atom_propagator(Real beta = 1.121,
			     Real mu0 = -.3){
    bool IsOk = true;

    
    ffd::lattice::HypercubicLattice<1> OneSite{{1}};

    
    ffd::imaginary_time::PeriodicImaginaryTime ThermalStrip(beta);
    
    
    auto O = CreateCoordinates(ThermalStrip,
			       OneSite);
    
    
    QuadraticAction S0 = Bar(Psi_(1)(O))*Psi_(1)(O)*(-mu0);
    S0 += FlipSpin(S0);
  

    
    auto G0 = ImaginaryTimeLatticePropagator(ThermalStrip,
    					     OneSite,
    					     S0);
    

    
    int N_samples = 100;
    for( int j: ffd::vector_range::Range(N_samples) ){
      

      Real tau_bpsi = 1e-15, tau_psi= (j+0.2)*beta/N_samples;

      
      auto psi = Psi_(-1)[0];
      auto Xpsi = O;
      component<0>(Xpsi) = tau_psi;
      std::any any_psi = Xpsi;
      auto psi_any = std::pair{psi, any_psi};


      auto bpsi = Bar(Psi_(-1))[0];
      auto Xbpsi = O;
      component<0>(Xbpsi) = tau_bpsi;
      std::any any_bpsi = Xbpsi;
      auto bpsi_any = std::pair{bpsi, any_bpsi};

    
      auto G0_exact = [mu0, beta](Real tau){return -exp(tau*mu0)/(1+exp(beta*mu0));};

    

      Real val = G0({{psi_any, bpsi_any}});

      
      // std::cerr<< val <<" "<< G0_exact(tau_psi-tau_bpsi) + val <<std::endl;
      IsOk = IsOk &&
	std::abs(G0_exact(tau_psi-tau_bpsi) + val) < 20*std::numeric_limits<Real>::epsilon();

      
    }

    
    return IsOk;
  }


}//namespace
