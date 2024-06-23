namespace ffd::imaginary_time_lattice_proposer::unit_test{

  void propose(){

    Real beta = 2.1231241;
    
    ffd::imaginary_time::PeriodicImaginaryTime beta_strip(beta);
    ffd::lattice::HypercubicLattice<2> chain(std::array<int, 2>{10, 12});
    auto O = ffd::user_space::CreateCoordinates(beta_strip, chain);
      
    
    auto gauss_one = [](Real tau){
	       return exp(-tau*tau/2);
	     };

    
    auto gauss_vec = [](auto x){
			 Real gauss_prod = 1;
			 for( int j: ffd::vector_range::Range( size(x) )){
			   gauss_prod *= std::exp(-.5*x[j]*x[j]);
			 }
			 return gauss_prod;
		       };

    
    Proposer proposer(beta_strip, chain,
		      gauss_one, gauss_vec);


    
    auto const N_samples = 1ul<<10;
    using ffd::vector_range::Range;
    for( auto j: Range(N_samples) ){
      auto time_prp = proposer.propose_time();
      // std::cerr<<j<<" "<<time_prp<<" ";
      

      assert( std::abs(time_prp) <= beta/2 );

      
      [[maybe_unused]] auto r = proposer.propose_space();
    //   for( int x: r ){
    // 	std::cerr<<x<<" ";
    //   }
      for( [[maybe_unused]] auto x: RealSpace( proposer.propose( O ) ) ) {
    	// std::cerr<<x<<" ";
      }

    //   std::cerr<<"\n";

    }
    



  }


}//namespace
