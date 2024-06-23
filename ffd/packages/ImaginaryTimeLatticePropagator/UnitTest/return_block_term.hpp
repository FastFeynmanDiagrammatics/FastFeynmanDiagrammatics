

namespace ffd::user_space::imaginary_time_lattice_propagator::unit_test{

  void return_block_term(){
    int const L = 10;
    ffd::lattice::HypercubicLattice<3> Cubic{{L, L, L}};

    

    Real const Beta = 2;
    ffd::imaginary_time::PeriodicImaginaryTime ThermalStrip(Beta);

    
    {
      auto O = CreateCoordinates(ThermalStrip, Cubic);
      auto X = O;
      component<1>(X) += 1;

    

      Real const t = 1;
      QuadraticAction S0 = Psi_(1)(O)*Bar(Psi_(1))(X)*(-t);
      S0 += HermitianConjugate(S0);
      Real const Mu = -1;
      S0 += Bar(Psi_(1))(O)*Psi_(1)(O)*(-Mu);
      S0 += FlipSpin(S0);


    
      ImaginaryTimeLatticePropagator G0{ThermalStrip,
					Cubic,
					S0};


    

      G0.ComputeInverseBravais();


    
      auto block_term0 =
	G0.ReturnBlockTerm(0);
      {
	auto [block_term_val, is_fermion, block_size, not_anomalous] = block_term0;

	static_assert((
		       std::is_same_v<decltype(block_term_val), std::vector<typename decltype(G0)::BlockTerm>>
		       ));
	
	// std::cerr<<"number elements = "<<std::size(block_term_val)<<std::endl;
	assert((
		std::size(block_term_val) == 3
		));

    
	for(int term: ffd::vector_range::Range(3)){
	  auto bilinear = block_term_val[term];
	  auto [value, R, indexes, daggers] = bilinear;


	  assert((
		  std::abs(value*value-1.) <
		  20*std::numeric_limits<Real>::epsilon()
		  ));


	  for(int dot: {0, 1}){
	    for(int j: ffd::vector_range::Range(3)){
	      assert((
		      std::abs( R[dot][j] - (dot == 1-term)*(j == 0)*(term != 2) ) <
		      20*std::numeric_limits<Real>::epsilon()
		      ));
	    }
	  }

      
	  for(int dot: {0, 1}){
	    assert((
		    indexes[dot] == 0
		    ));
	  }
      

	  for(int dot: {0, 1}){
	    assert((
		    daggers[dot] == (term < 2 ? dot : 1-dot)
		    ));
	  }
      
	}
    
	// std::cerr<<"is_fermion = "<<std::boolalpha<<is_fermion<<std::endl;
	assert((
		is_fermion
		));
	// std::cerr<<"block_size = "<<block_size<<std::endl;
	assert((
		block_size == 1
		));
	// std::cerr<<"not_anomalous has value = "<<not_anomalous.has_value()<<std::endl;
	assert((
		!not_anomalous.has_value()
		));
      }


    
      auto block_term1 =
	G0.ReturnBlockTerm(1);

      {
	auto [block_term_val, is_fermion, block_size, not_anomalous] = block_term1;

	// std::cerr<<"number elements = "<<std::size(block_term_val)<<std::endl;
	assert((
		std::size(block_term_val) == 3
		));

    
	for(int term: ffd::vector_range::Range(3)){
	  auto bilinear = block_term_val[term];
	  auto [value, R, indexes, daggers] = bilinear;


	  assert((
		  std::abs(value*value-1.) <
		  20*std::numeric_limits<Real>::epsilon()
		  ));


	  for(int dot: {0, 1}){
	    for(int j: ffd::vector_range::Range(3)){
	      assert((
		      std::abs( R[dot][j] - (dot == 1-term)*(j == 0)*(term != 2) ) <
		      20*std::numeric_limits<Real>::epsilon()
		      ));
	    }
	  }

      
	  for(int dot: {0, 1}){
	    assert((
		    indexes[dot] == 0
		    ));
	  }
      

	  for(int dot: {0, 1}){
	    assert((
		    daggers[dot] == (term < 2 ? dot : 1-dot)
		    ));
	  }
      
	}
    
	// std::cerr<<"is_fermion = "<<std::boolalpha<<is_fermion<<std::endl;
	assert((
		is_fermion
		));
	// std::cerr<<"block_size = "<<block_size<<std::endl;
	assert((
		block_size == 1
		));
	// std::cerr<<"not_anomalous has value = "<<not_anomalous.has_value()<<std::endl;
	assert((
		!not_anomalous.has_value()
		));
      }
    }


    




    ffd::lattice::KagomeLattice K{{L+1, L+3}};

    {

      constexpr int d = 2;
      

      auto O = CreateCoordinates(ThermalStrip, K);

      
      Complex Mu = -2;
      QuadraticAction S0 = Bar(Eta_(1)(O))*Eta_(1)(O)*(-Mu);

      
      auto X = O;
      Complex t = 1. + ffd::user_space::I;


      component<d+1>(X) = 1;
      S0 += Bar(Eta_(1)(O))*Eta_(1)(X)*(-t);


      component<d+1>(X) = 2;
      S0 += Bar(Eta_(1)(O))*Eta_(1)(X)*(-t);


      component<d>(X) += 1;
      S0 += Bar(Eta_(1)(O))*Eta_(1)(X)*(-t);
      
      component<d-1>(X) += 1;
      S0 += Bar(Eta_(1)(O))*Eta_(1)(X)*(-t);



      S0 += HermitianConjugate(S0);


      Complex lambda = -1. + 2.*ffd::user_space::I;
      S0 += Bar(Eta_(1)(O))*Eta_(-1)(X)*lambda;

      
      S0 += FlipSpin(S0);



      ImaginaryTimeLatticePropagator G0{ThermalStrip,
					K,
					S0};


      
      G0.ComputeInverseBravais();
      

      auto block0 = G0.ReturnBlockTerm(0);


      
      auto [block_terms,
	    is_fermion,
	    block_size,
	    not_anomalous] = block0;



      // std::cerr<<"number_of_terms = "<<std::size(block_terms)<<std::endl;
      assert((
	      std::size(block_terms) == 22
	      ));

      // std::cerr<<"is_fermion = "<<std::boolalpha<<is_fermion<<std::endl;
      assert((
	      !is_fermion
	      ));
      // std::cerr<<"block_size = "<<block_size<<std::endl;

      assert((
	      block_size == 6
	      ));

    }



    ffd::lattice::HypercubicLattice<1> Chain{{L}};


    {

      
      auto O = CreateCoordinates(ThermalStrip,
				 Chain);

      
      auto X = O;
      component<1>(X) += 1;

      Complex t = ffd::user_space::I;
      
      QuadraticAction S0 = Rho_(0)(O)*Rho_(0)(X)*t;



      ImaginaryTimeLatticePropagator G0{ThermalStrip,
					Chain,
					S0};

      
      G0.ComputeInverseBravais();


      auto block0 = G0.ReturnBlockTerm(0);

      
      auto [block_terms,
	    is_fermion,
	    block_size,
	    not_anomalous] = block0;


      // std::cerr<<"number of terms = "<<std::size(block_terms)<<std::endl;
      assert((
	      std::size(block_terms) == 4
	      ));

      // std::cerr<<"is_fermion = "<<is_fermion<<std::endl;
      assert((
	      is_fermion
	      ));
      // std::cerr<<"block_size = "<<block_size<<std::endl;
      assert((
	      block_size == 1
	      ));

      for(int j: ffd::vector_range::Range(4)){
	assert((
		std::abs(std::get<0>(block_terms[j]) - ffd::user_space::I/2.) <
		20*std::numeric_limits<Real>::epsilon()
		));
      }

    }
    
  }

}//namespace
