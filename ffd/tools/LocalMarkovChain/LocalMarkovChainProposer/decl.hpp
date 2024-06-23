

namespace ffd::local_markov_chain{

  struct VoidStruct{};

  
  
  template<bool IndependentProposals,
	   typename ProposerType>
  
  struct LocalMarkovChainProposer{

    
    ProposerType Proposer;

    
    LocalMarkovChainProposer(ProposerType proposer_):
      Proposer(proposer_) {}

    LocalMarkovChainProposer(ProposerType proposer_,
			     bool):
      Proposer(proposer_) { }

    
    using PositionVariable = int;
    using WeightOfMove = Real;

    
    template<typename VariableType,
	     typename ConstantType = VoidStruct>
    
    std::pair<std::pair<PositionVariable,
			VariableType>,
	      WeightOfMove>
    ProposeMove(std::vector<VariableType> variables_,
		std::vector<ConstantType> constants_ =
		std::vector<ConstantType>(),
		int which_one_to_move = -1);


    
    template<int num_moves,
	     typename VariableType,
	     typename ConstantType = VoidStruct>

    std::pair<std::array<std::pair<PositionVariable,
				   VariableType>,
			 num_moves>,
	      WeightOfMove>
    ProposeMoves(std::vector<VariableType> variables_,
		 std::vector<ConstantType> constants_ =
		 std::vector<ConstantType>());

  };

  
  template<typename T>
  LocalMarkovChainProposer(T x) -> LocalMarkovChainProposer<true, T>;


  template<typename T>
  LocalMarkovChainProposer(T x, bool) -> LocalMarkovChainProposer<false, T>;

	
}//namespace
