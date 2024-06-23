

namespace ffd::local_markov_chain{

  
  template<bool IndependentProposals,
	   typename ProposerType>
  template<typename FuncReturningTwoValues_t,
	   typename VariableType,
	   typename ConstantType,
	   typename FuncReturnType>

  
  std::tuple<bool,
	     std::array<FuncReturnType, 2>,
	     Real,
	     std::pair<int, VariableType>>

  
  LocalMarkovChain<IndependentProposals,
		   ProposerType>::
  Move(std::array<FuncReturnType, 2> func_now,
       FuncReturningTwoValues_t const& Func,
       std::vector<VariableType> variables,
       std::vector<ConstantType> constants){

    
    if( NeverMoved ){
      if constexpr (std::is_same_v<ConstantType, VoidStruct>){
	func_now = Func(variables);
      }else{
	func_now = Func(variables, constants);
      }
      NeverMoved = false;
      auto proba_first = MixingProbabilityDistributions(func_now, LambdaMixing);
      return std::tuple{true, func_now, proba_first, std::pair{0, variables[0]}};
    }
    
    
    
    auto [position_value, weight] =
      Proposer.ProposeMove(variables,
			   constants);

    
    
    auto [position, value] = position_value;
    variables[position] = value;


    std::array<FuncReturnType, 2> func_prp;
    if constexpr(std::is_same_v<ConstantType, VoidStruct>){
	func_prp = Func(variables);
      }else{
      func_prp = Func(variables, constants);
    }
		  

    


    Real proba_now = ffd::local_markov_chain::
      MixingProbabilityDistributions(func_now, LambdaMixing);
    Real proba_prp = ffd::local_markov_chain::
      MixingProbabilityDistributions(func_prp, LambdaMixing);


    bool IsAccepted = ffd::random_distributions::
      Proba() < proba_prp*weight/proba_now;


    return std::tuple{IsAccepted, func_prp, proba_prp, position_value};
  }

}//namespace
