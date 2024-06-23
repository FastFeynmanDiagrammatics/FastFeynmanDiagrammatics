

namespace ffd::eigenvalues_vectors{

  template<typename Field>
  auto
  MultiplyMatrices(std::vector<std::vector<Field>> const& factor1_,
		   std::vector<std::vector<Field>> const& factor2_){


    using std::size;
    

    assert( size(factor1_) == size(factor2_) );
    auto const size_matrices = size(factor1_);
    std::vector<std::vector<Field>>
      result(size_matrices, std::vector<Field>(size_matrices, 0.));

    
    for(std::size_t j=0; j < size_matrices; ++j){
      assert( size(factor1_[j]) == size_matrices);
      assert( size(factor2_[j]) == size_matrices);
      for(std::size_t k=0; k < size_matrices; ++k){
	for(std::size_t m=0; m < size_matrices; ++m){
	  result[j][k] += factor1_[j][m]*factor2_[m][k];
	}
      }
    }
    return result;
  }

}//namespace
