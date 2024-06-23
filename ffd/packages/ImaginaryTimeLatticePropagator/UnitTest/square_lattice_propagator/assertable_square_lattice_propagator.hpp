

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  auto
  assertable_square_lattice_propagator(int Lx = 8,
				       int Ly = 7,
				       Real tx = 1.21312,
				       Real ty = .872121,
				       Real txy = .2141231,
				       Real mu0 = -.812312,
				       Real beta = 1.812){
    bool IsOk = true;
    

    
    ffd::lattice::HypercubicLattice<2> SquareLattice{{Lx, Ly}};
    ffd::imaginary_time::PeriodicImaginaryTime BetaStrip(beta);


    
    auto O = CreateCoordinates(BetaStrip,
			       SquareLattice);


    
    QuadraticAction S0 = Bar(Psi_(1)(O))*Psi_(1)(O)*(-mu0);
    auto X = O;

    
    component<1>(X) = 1;
    auto hop = Bar(Psi_(1)(O))*Psi_(1)(X)*(-tx);
    S0 += hop + HermitianConjugate(hop);

    
    component<2>(X) = 1;
    hop = Bar(Psi_(1)(O))*Psi_(1)(X)*(-txy);
    S0 += hop + HermitianConjugate(hop);


    component<2>(X) = -1;
    hop = Bar(Psi_(1)(O))*Psi_(1)(X)*(-txy);
    S0 += hop + HermitianConjugate(hop);


    component<1>(X) = 0;
    component<2>(X) = 1;
    hop = Bar(Psi_(1)(O))*Psi_(1)(X)*(-ty);
    S0 += hop + HermitianConjugate(hop);


    
    S0 += FlipSpin(S0);


    
    auto G0 = ImaginaryTimeLatticePropagator(BetaStrip,
					     SquareLattice,
					     S0);

    

    
    auto G0_exact = [mu0, beta, tx, ty, txy, Lx, Ly](int x, int y, Real tau){
		      Real ret = 0;
		      for(int nx=0; nx < Lx; ++nx){
			
			Real kx = 2*nx*ffd::core_math::Pi/Lx;
			
			for(int ny=0; ny < Ly; ++ny){
			  Real ky = 2*ny*ffd::core_math::Pi/Ly;
			  
			  Real eps_k = -2*tx*cos(kx)-2*ty*cos(ky)-4*txy*cos(kx)*cos(ky)-mu0;
			  
			  ret += -cos(ky*y)*cos(kx*x)*exp(-tau*eps_k)/(1+exp(-beta*eps_k));
			  
			}
		      }
		      return ret/Lx/Ly;
		    };



    auto psi = Psi_(-1)[0], bpsi = Bar(Psi_(-1))[0];
    Real tau0 = 1e-15;
    component<0>(O) = tau0;
    auto field_pos1 = std::pair{bpsi, std::any(O)};
    

    
    Real N_samples = 30; 
    for(int x=0; x < Lx; ++x){
      for(int y=0; y < Ly; ++y){
	for(int j=0; j<N_samples; ++j){
	  Real tau = (j+.2)*beta/N_samples+tau0;

	  component<0>(X) = tau;
	  component<1>(X) = x;
	  component<2>(X) = y;

	  auto field_pos0 = std::pair{psi, std::any(X)};
	  // std::cerr<<G0_exact(x, y, tau)<<" "<<G0_exact(x, y, tau)+G0({{field_pos0, field_pos1}})<<std::endl;
	  IsOk = IsOk &&
	    std::abs(G0_exact(x, y, tau)+G0({{field_pos0, field_pos1}})) <
	    20*std::numeric_limits<Real>::epsilon();
	}
      }
    }
        
    

    return IsOk;
  }
  



}//namespace
