

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void test_hamiltonian_k(){
    Real precision = 1e-12;

    
    int Lx = 13, Ly=11;
    ffd::lattice::HypercubicLattice<2> Square{{Lx, Ly}};

    Real Beta = 3.;
    ffd::imaginary_time::PeriodicImaginaryTime ImagTime(Beta);

    
    auto O = CreateCoordinates(ImagTime, Square);
    auto Y = O;

    Real Mu = -2.28;
    QuadraticAction Action0 = Bar(Psi_(1)(O))*Psi_(1)(Y)*(-Mu);
    
    
    component<1>(Y) += 1;
    Real t_x = 2.123123;
    auto hopping_x = Bar(Psi_(1)(Y))*Psi_(1)(O)*(-t_x);
    Action0 += hopping_x + HermitianConjugate(hopping_x);
    
    
    
    Y = O;
    component<2>(Y) += 1;
    Real t_y = 1.253798;
    auto hopping_y = Bar(Psi_(1)(Y))*Psi_(1)(O)*(-t_y);
    Action0 += hopping_y + HermitianConjugate(hopping_y);


    
    auto G0 = ImaginaryTimeLatticePropagator(ImagTime,
  					     Square,
  					     Action0);

        
    

    auto [block_terms, is_fermion, block_size, not_anomalous]
      = G0.ReturnBlockTerm(0);

    
    auto eps_k = [Mu, t_x, t_y](Real kx, Real ky)
  		 {return -Mu -2*t_x*std::cos(kx)-2*t_y*std::cos(ky);};

    
    
    for(int j=0; j < Lx; ++j){
      Real k_x = 2*j*ffd::core_math::Pi/Lx;
      for(int k=0; k < Ly; ++k){
  	Real k_y = 2*k*ffd::core_math::Pi/Ly;
  	auto H = G0.ReturnBlockHamiltonian_k(block_terms,
					     is_fermion,
					     block_size,
					     {{k_x, k_y}});
	assert((
		std::size(H) == 2
		));
	for(int u: {0, 1}){
	  for(int v: {0, 1}){
	    //   std::cerr<<std::setprecision(15)<<
	    // "("<<k_x<<", "<<k_y<<") = "<<
	    // eps_k(k_x, k_y)<<" "<<H[u][v]<<std::endl;

	  assert((
	  	  std::abs( (1-2*u)*(u==v)*eps_k(k_x, k_y) - H[u][v] ) <
	  	  precision
  	  	));
	  
	  }
	}
      }
    }

  }

}//namespace
