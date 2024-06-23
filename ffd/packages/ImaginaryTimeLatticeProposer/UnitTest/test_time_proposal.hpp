#pragma once

namespace ffd::imaginary_time_lattice_proposer::unit_test{

  void test_time_proposal(){

    Real beta = 2.1231241;
    
    ffd::imaginary_time::PeriodicImaginaryTime beta_strip(beta);
    
    ffd::lattice::HypercubicLattice<1> chain(std::array<int, 1>{10});

    auto f = [](Real tau){
	       return exp(-tau*tau);
	     };

    auto void_lambda = [](auto){return 1.; };

    Proposer proposer(beta_strip, chain,
		      f, void_lambda);

    
    auto const N_samples = 1ul<<10;
    using ffd::vector_range::Range;
    for( auto j: Range(N_samples) ){
      auto time_prp = proposer.propose_time();
      // std::cerr<<j<<" "<<time_prp<<std::endl;
      assert( std::abs(time_prp) <= beta/2 );
    }
    

    // auto r = proposer.propose_space();
    // for( int x: r ){
    //   std::cerr<<x<<" ";
    // }
    
  }

}//namespace
