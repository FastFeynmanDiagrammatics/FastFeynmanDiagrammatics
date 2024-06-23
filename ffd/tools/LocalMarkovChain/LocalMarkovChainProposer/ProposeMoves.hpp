

namespace ffd::local_markov_chain{

  template<bool IndependentProposals,
	   typename ProposerType>
  template<int n_moves,
	   typename VariableType,
	   typename ConstantType>

  std::pair<std::array<std::pair<int,
				 VariableType>,
		       n_moves>,
	    Real>
  
  LocalMarkovChainProposer<IndependentProposals,
			   ProposerType>::
  ProposeMoves(std::vector<VariableType> variables_,
	       std::vector<ConstantType> constants_){
    Real weight = 1;

    
    std::array<std::pair<int, VariableType>,
	       n_moves> PairsVariableMoved;

    
    std::set<int> variables_moved;

    
    for( int u: ffd::vector_range::Range(n_moves) ){

      
      int variable_to_move;
      do{
	variable_to_move = ffd::random_distributions::
	  RandomIntFromZeroTo( size(variables_)-1 );
      }while(  variables_moved.count( variable_to_move )  ==  1  );
      variables_moved.insert(variable_to_move);


      auto [pair_moved, weight_move] = ProposeMove(variables_,
						   constants_,
						   variable_to_move);
      auto [position, value] = pair_moved;
      variables_[position] = value;
      
      
      weight *= weight_move;
      PairsVariableMoved[u] = pair_moved;
      
    }
    
    return std::pair{PairsVariableMoved, weight};
  }

}//namespace
