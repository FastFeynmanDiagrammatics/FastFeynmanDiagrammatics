

namespace ffd::eigenvalues_vectors{

  void
  Hermitize(std::vector<std::vector<Complex>>& H_){


    int const size = std::size(H_);
    auto H_copy = H_;


    for(int j=0; j < size; ++j){
      for(int k=0; k < size; ++k){
	H_[j][k] = .5*(H_copy[j][k] + std::conj(H_copy[k][j]));
      }
    }
  }
  


}//namespace
