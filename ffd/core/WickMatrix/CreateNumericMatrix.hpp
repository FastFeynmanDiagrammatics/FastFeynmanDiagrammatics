namespace ffd::user_space{

  template<typename Field>
  
  auto
  CreateNumericMatrix(ffd::wick_matrix::
		      WickMatrix<ffd::feynman_edge::FeynmanEdge>
		      const& EdgeMatrix,
		      ffd::feynman_edge::FeynmanEdgeMap<Field> const& Map_){
    std::size_t const SizeEdgeMatrix = size( EdgeMatrix.Components );

    
    std::vector<Field> ret(SizeEdgeMatrix);

    
    for(std::size_t k = 0ul; k < SizeEdgeMatrix; ++k){
      ret[k] = Map_.at( EdgeMatrix.Components[k] );//Map_[EdgeMatrix.Components[k]];
    }

    
    return ret;
  }



}//namespace
