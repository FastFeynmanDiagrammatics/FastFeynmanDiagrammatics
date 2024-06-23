namespace ffd::user_space{
  
  template<int order, typename coord_t>
  std::array<std::array<int, (1<<order)>, 2>
  compute_signs_from_chis(std::vector<QuantumFieldPair<coord_t>> const& chis){
    std::array<std::array<int, (1<<order)>, 2> ret;
    for(int j: {0, 1}){
      ret[j].fill(1);
    }
    
    
    for(uint S=0; S<(1<<order); ++S){
      auto vec_S = ffd::set_theory::
	VectorOfBinaryDigitsOf(S);
      for(uint j=0; j<size(vec_S); ++j){
	for(int s: {0, 1}){
	  ret[s][S] *= chis[vec_S[j]].sign;
	}
      }
      for(uint j=2*order; j<size(chis); ++j){
	ret[1][S] *= chis[j].sign;
      }
    }
    
    
    return ret;
  }

}//namespace
