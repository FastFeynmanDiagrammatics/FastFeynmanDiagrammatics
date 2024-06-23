namespace ffd::user_space{
  
  template<typename intType = ffd::wick_matrix::intTypeDefault>
  
  std::vector<ffd::wick_matrix::WickMatrix<FeynmanEdge, intType>>
  
  CreateFeynmanEdgeWickMatrices(ffd::wick_function::WickFunction const& W_){
    std::vector<ffd::wick_matrix::WickMatrix<FeynmanEdge, intType>> ret(size(W_));

    
    for(std::size_t j=0; j < size(W_); ++j){
      auto [Block0, Block1] = W_[j];//can optimize by using pointers
      assert(size(Block0) != 0);
      assert(size(Block1) == 0 || size(Block1) == size(Block0));
      if(size(Block1) == 0){
	Block1 = Block0;//can optimize by using pointers
      }

      
      auto Qfield = Block0[0].first;
      ret[j] = ffd::wick_matrix::WickMatrix<FeynmanEdge>(Block0[0].first, size(Block0));
      for(std::size_t m=0; m < size(Block0); ++m){
	for(std::size_t k=(Qfield.Dagger() == 0)*(m+Qfield.IsFermion()); k < size(Block0); ++k){
	  auto QPm = Block0[m];
	  auto QPk = Block1[k];
	  auto QPmin = QPm < QPk ? QPm : QPk;
	  auto QPMax = QPm < QPk ? QPk : QPm;
	  ret[j](m, k, 'a') = {QPmin, QPMax};
	}
      }
    }
    return ret;
  }
  
}//namespace ffd::user_space


