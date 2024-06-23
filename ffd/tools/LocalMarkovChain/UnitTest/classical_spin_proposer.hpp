

namespace ffd::local_markov_chain::unit_test{

  struct func_test{
    std::array<Real, 2>
    operator()(std::vector<bool>) const{
      return std::array{1., 1.};
    }
  };

  
  void classical_spin_proposer(){

    ClassicalSpinProposer spin_proposer;

    std::vector<bool> spins{true, true, true, true};
    
    
    LocalMarkovChainProposer Chain{spin_proposer};

    
    LocalMarkovChain Markov{Chain, .1};

    func_test Func;
    
    for( [[maybe_unused]] int j: ffd::vector_range::Range(10) ){
      [[maybe_unused]] auto [pair, weight] =
	Chain.ProposeMove(spins);

      [[maybe_unused]] auto [is_accepted, dist, proba_prp, pos_value] = 
	Markov.Move({0., 0.}, Func, spins);

      [[maybe_unused]] auto [spin_number, value] = pair;

      spins[spin_number] = value;
      
      // std::cerr<<spin_number<<" "<<value<<" "<<weight<<std::endl;
    }

    // std::cerr<<"\n";
    for( [[maybe_unused]] int j: ffd::vector_range::Range(10) ){
      [[maybe_unused]]
	auto [pair_array, weight] =
	Chain.ProposeMoves<2>(spins);

      for( int j: {0, 1} ){
	auto [spin_number, value] = pair_array[j];
	spins[spin_number] = value;
	// std::cerr<<spin_number<<" "<<value<<" "<<weight<<std::endl;
      }

    }
    
  }

}//namespace
