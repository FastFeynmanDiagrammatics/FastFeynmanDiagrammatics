

namespace ffd::eigenvalues_vectors{

  auto
  Dagger(std::vector<std::vector<Complex>> const& M_){
    int const size = std::size(M_);
    std::vector<std::vector<Complex>> M_dagger(size, std::vector<Complex>(size, 0.));


    
    for(int j = 0; j < size; ++j){
      for(int k = 0; k < size; ++k){
	M_dagger[j][k] = std::conj(M_[k][j]);
      }
    }
    return M_dagger;
  }


}//namespace
