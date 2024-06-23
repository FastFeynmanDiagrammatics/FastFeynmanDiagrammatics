namespace ffd::wick_function{

  WickFunction& WickFunction::operator*=(std::pair<QuantumFieldVertex, char> QV_){
    using std::size;
    auto [V_, j_] = QV_;
    for(char dot=0; dot < char(size(V_)); ++dot){
      std::array<char, 2> VertexDot{j_, dot};
      const auto& Dot = V_[dot];
      for(auto const& QF: Dot){
	QuantumFieldVertexDot QP_{QF, VertexDot};
	const auto WhatIsMyBlock = WhichBlock(QF);
	const auto BlockNumber_it = std::find(BlockNumbers.begin(), BlockNumbers.end(), WhatIsMyBlock);
	std::size_t BlockNumber = std::distance(BlockNumbers.begin(), BlockNumber_it);
	if( BlockNumber != size(*this) ){
	  if(QF.IsFermion()){
	    int NumberPassedFermions = 0;
	    for(std::size_t j = size(*this)-1; j > BlockNumber - (QF.Dagger() == -1); --j){
	      const auto& [Block0_1, Block1] = (*this)[j];
	      if(IsFermionicBlock[j]){
		NumberPassedFermions += size(Block0_1);
		NumberPassedFermions += size(Block1);
	      }
	    }
	    Sign *= 1 - 2*(NumberPassedFermions%2);
	  }
	}else{
	  (*this).resize(BlockNumber+1);
	  if(QF.Dagger() != 0){
	    (*this)[BlockNumber][1] = std::vector<QuantumFieldVertexDot>();
	  }
	  BlockNumbers.resize(BlockNumber+1);
	  BlockNumbers[BlockNumber] = WhatIsMyBlock;
	  IsFermionicBlock.resize(BlockNumber+1);
	  IsFermionicBlock[BlockNumber] = QF.IsFermion();
	}
	auto& [Block0_1, Block1] = (*this)[BlockNumber];
	if(QF.Dagger() <= 0){
	  Block0_1.push_back(QP_);
	}else{
	  Block1.push_back(QP_);
	}
      }
    }
    return *this;
  }

}//namespace ffd::wick_function
