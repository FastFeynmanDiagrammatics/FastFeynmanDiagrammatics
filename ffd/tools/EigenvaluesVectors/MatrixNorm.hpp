

namespace ffd::eigenvalues_vectors{

  template<int n, typename Field>
  auto MatrixNorm(std::vector<std::vector<Field>> M_){
    Real norm_n = 0.;

    
    
    for(std::size_t j = 0; j < size(M_); ++j){
      for(std::size_t k = 0; k < size(M_); ++k){
	norm_n += std::pow(std::abs(M_[j][k]), n);
      }
    }
    return std::pow(norm_n, 1./n);
  }


}//namespace
