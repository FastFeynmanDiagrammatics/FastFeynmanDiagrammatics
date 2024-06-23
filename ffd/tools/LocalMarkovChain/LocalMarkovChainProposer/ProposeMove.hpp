

namespace ffd::local_markov_chain{

  template<bool IndependentProposals,
	   typename ProposerType>

  template<typename VariableType,
	   typename ConstantType>

  std::pair<std::pair<int, VariableType>,
	    Real>

  LocalMarkovChainProposer<IndependentProposals,
			   ProposerType>::
  ProposeMove(std::vector<VariableType> variables_,
	      std::vector<ConstantType> constants_,
	      int which_one){
    int variable_number_to_move = which_one;

    
    if(which_one < 0){
      variable_number_to_move = ffd::random_distributions::
      RandomIntFromZeroTo( size(variables_)-1 );
    }

    
    std::pair<VariableType, Real> return_proposer;

    
    if constexpr ( IndependentProposals ){
	VariableType return_independent = Proposer(variables_[variable_number_to_move]);
	return_proposer = std::pair{return_independent, Real(1.)};
      }else if constexpr ( std::is_same_v<ConstantType, VoidStruct> ){
return_proposer = Proposer(variable_number_to_move,
			   variables_);
}else{
      return_proposer =  Proposer(variable_number_to_move,
				  variables_,
				  constants_);
    }

    
    auto [variable, weight] = return_proposer;
    return std::pair{std::pair{variable_number_to_move, variable}, weight};
  }

  
}//namespace
