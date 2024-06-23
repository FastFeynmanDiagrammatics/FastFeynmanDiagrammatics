namespace ffd::user_space{

  template<typename Field,
	   typename intType = ffd::wick_matrix::intTypeDefault>
  
  ffd::wick_matrix::WickMatrix<Field>
  
  CreateNumericWickMatrix(ffd::wick_matrix::WickMatrix<FeynmanEdge> const& EdgeMatrix,
			  ffd::feynman_edge::FeynmanEdgeMap<Field> const& Map_){
    ffd::wick_matrix::WickMatrix<Field> ret;
    const intType SizeEdgeMatrix = size(EdgeMatrix.Components);
    ret.Components.resize(SizeEdgeMatrix);
    ret.NumberDestructionOperators = EdgeMatrix.NumberDestructionOperators;
    // ret.IsFermion = EdgeMatrix.IsFermion;
    // ret.NotHermitian = EdgeMatrix.NotHermitian;
    for(intType k=0; k < SizeEdgeMatrix; ++k){
      ret.Components[k] = Map_.at(EdgeMatrix.Components[k]);//Map_[EdgeMatrix.Components[k]];
    }
    return ret;
  }

}//namespace ffd::user_space
