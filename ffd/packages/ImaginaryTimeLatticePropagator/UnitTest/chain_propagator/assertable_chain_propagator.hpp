

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  auto
  assertable_chain_propagator(int L = 11,
			      Real t = 1.321,
			      Real mu0 = -.612312,
			      Real beta = 2.12){
    bool IsOk = true;
    

    
    ffd::lattice::HypercubicLattice<1> Chain{{L}};
    ffd::imaginary_time::PeriodicImaginaryTime BetaStrip(beta);


    
    auto O = CreateCoordinates(BetaStrip,
			       Chain);


    
    QuadraticAction S0 = Bar(Psi_(1)(O))*Psi_(1)(O)*(-mu0);
    auto X = O;
    component<1>(X) += 1;
    auto hop = Bar(Psi_(1)(O))*Psi_(1)(X)*(-t);
    S0 += hop + HermitianConjugate(hop);
    S0 += FlipSpin(S0);


    
    auto G0 = ImaginaryTimeLatticePropagator(BetaStrip,
					     Chain,
					     S0);

    

    
    auto G0_exact = [mu0, beta, t, L](int x, Real tau){
		      Real ret = 0;
		      for(int n=0; n < L; ++n){
			Real k = 2*n*ffd::core_math::Pi/L;
			Real eps_k = -2*t*cos(k)-mu0;
			ret += -cos(k*x)*exp(-tau*eps_k)/(1+exp(-beta*eps_k));
		      }
		      return ret/L;
		    };



    auto psi = Psi_(-1)[0], bpsi = Bar(Psi_(-1))[0];
    Real tau0 = 1e-15;
    component<0>(O) = tau0;
    auto field_pos1 = std::pair{bpsi, std::any(O)};
    

    
    Real N_samples = 100; 
    for(int x=0; x < L; ++x){
      for(int j=0; j<N_samples; ++j){
	Real tau = (j+.2)*beta/N_samples+tau0;

	component<0>(X) = tau;
	component<1>(X) = x;

	auto field_pos0 = std::pair{psi, std::any(X)};
	// std::cerr<<G0_exact(x, tau)<<" "<<G0_exact(x, tau)+G0({{field_pos0, field_pos1}})<<std::endl;
	IsOk = IsOk &&
	  std::abs(G0_exact(x, tau)+G0({{field_pos0, field_pos1}})) <
	  20*std::numeric_limits<Real>::epsilon();
      }
    }
        
    

    return IsOk;
  }
  


}//namespace
