namespace ffd::user_space{
  
  template<typename G0_t,
	   typename coord_t>
  auto
  fill_G0_matrix(G0_t const& G0,
		 std::vector<QuantumFieldPair<coord_t>> chis){
    std::vector<decltype(G0(chis[0]))> matrix(size(chis)*size(chis));
    for(uint j=0; j<size(chis); ++j){
      for(uint u=0; u<size(chis); ++u){
	matrix[u+j*size(chis)] =
	  G0(QuantumFieldPair<coord_t>(chis[j].fields[0], chis[u].fields[1]));
      }
    }
    
    
    return matrix;
  }

}//namespace
