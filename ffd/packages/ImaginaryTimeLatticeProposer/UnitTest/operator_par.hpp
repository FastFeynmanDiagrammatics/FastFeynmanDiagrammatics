namespace ffd::imaginary_time_lattice_proposer::unit_test{

  void operator_par(){
    using ffd::vector_range::Range;
    
    Real beta = 28;

    
    ffd::imaginary_time::PeriodicImaginaryTime beta_strip(beta);
    int Lx = 11, Ly=12;
    ffd::lattice::HoneycombLattice H({Lx, Ly});
    auto O = ffd::user_space::CreateCoordinates(beta_strip, H);

    

    Real sigma_time = 1., sigma_space = 1;

    
    
    auto gauss_one = [=](Real tau){
		       return std::exp(-.5*std::pow(tau/sigma_time, 2));
	     };

    
    
    auto gauss_vec = [=](auto x){
			 Real gauss_prod = 1;
			 for( int j: ffd::vector_range::Range( size(x) ) ){
			   gauss_prod *= std::exp(-.5*std::pow(x[j]/sigma_space, 2));
			 }
			 return gauss_prod;
		       };

    
    
    
    Proposer proposer(beta_strip, H,
		      gauss_one, gauss_vec);

    
    
    auto [coordinates,
    	  int_components] = ffd::lattice::Range(H);
    assert(( size(int_components) == size(coordinates) ));
    assert(( size(int_components) == Lx*Ly*2 ));


    assert((
	    std::abs( proposer.time_cumulative(-beta/2) ) < 1e-10
	    ));
    assert((
	    std::abs( proposer.time_cumulative(beta/2) - 1) < 1e-10
	    ));


    constexpr int n_atoms = decltype(H)::number_atoms_unit_cell; 
    for( auto atom: Range( n_atoms ) ){
      assert((
	      std::abs( proposer.space_cumulative[atom][size(coordinates)-1] -1) < 1e-10
	      ));
    }

    


    auto const N_samples_tau = 1ul<<5;
    for( auto j: Range(N_samples_tau) ){
      Real tau0 = -beta/2 + (j+.2)*beta/N_samples_tau;
      auto X = O;
      ffd::user_space::component<0>(X) += tau0;
      // std::cerr<<tau0<<" "<<proposer(X,O)[0]/proposer(O,O)[0]-
      // 	gauss_one(tau0)<<std::endl;

      
      assert((
      	      std::abs( proposer(X,O)/proposer(O,O)-
      			gauss_one(tau0) ) < 1e-10
      	      ));
    }


    

    
    // for( auto atom: Range(n_atoms) ){
    //   for(  auto  j:  Range( size(int_components) )  ){
    // 	auto [x, y, z] = int_components[j];
    // 	std::cerr<<atom<<": "<<x<<" "<<y<<" "<<z<<" ";
    // 	auto [x1, y1] = RealSpace(coordinates[j]);
    // 	std::cerr<<x1<<" "<<y1<<" "<<proposer.space_function[atom][j]<<
    // 	  " "<<proposer.space_cumulative[atom][j]<<"\n";
    //   }
    // }


    
	  
    
    auto const N_samples = 1ul<<8;
    for( auto j: Range(N_samples) ){
      auto time_prp = proposer.propose_time();
      // std::cerr<<j<<" "<<time_prp<<" ";
      

      assert( std::abs(time_prp) <= beta/2 );

      
      [[maybe_unused]] auto r = proposer.propose_space();
      // for( int x: r ){
      // 	std::cerr<<x<<" ";
      // }
      // std::cerr<<",  X = ";
      auto X = proposer.propose( O );
      // std::cerr<<ffd::get<0>(X)()<<" ";
      for(  [[maybe_unused]]  auto  x:  RealSpace( X )  ) {
    	// std::cerr<<x<<" ";
      }
      // std::cerr<<", P = ";
      auto xo = proposer(X, O);
      auto ox = proposer(O, X);
      // std::cerr<<xo<<" "<<ox<<" ";
      
      // std::cerr<<"\n";
      assert((
	      std::abs( xo - ox ) < 1e-10
	      ));
      
    }

    


  }


}//namespace
