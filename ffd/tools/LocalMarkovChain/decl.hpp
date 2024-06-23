

namespace ffd::local_markov_chain{

  template<bool IndependentProposals,
	   typename ProposerType>
  
  struct LocalMarkovChain{

    LocalMarkovChainProposer<
      IndependentProposals,
      ProposerType> Proposer;

    Real LambdaMixing;

    bool NeverMoved;
    
    
    LocalMarkovChain(LocalMarkovChainProposer<
		     IndependentProposals,
		     ProposerType> Proposer_,
		     Real LambdaMixing_,
		     bool NeverMoved_ = true):
      Proposer(Proposer_),
      LambdaMixing(LambdaMixing_),
      NeverMoved(NeverMoved_) {}
    
    
    

    using IsAccepted_t = bool;
    using VariablePosition = int;
    using probability_t = Real;
    

    template<typename FuncReturningTwoValues_t,
	     typename VariableType,
	     typename ConstantType = VoidStruct,
	     typename FuncReturnType = Real>
    
    std::tuple<IsAccepted_t,
	       std::array<FuncReturnType, 2>,
	       probability_t,
	       std::pair<VariablePosition,
			 VariableType>>
    
    Move(std::array<FuncReturnType, 2> func_now,
	 FuncReturningTwoValues_t const& Func,
	 std::vector<VariableType> variables_now,
	 std::vector<ConstantType> constants_now =
	 std::vector<ConstantType>());
    
    
  };

}//namespace
