

namespace ffd::eigenvalues_vectors{

  template<int n, typename Field>
  auto TraceNorm(std::vector<std::vector<Field>> M_){
    Real norm_n = 0.;

    
    for(std::size_t j = 0; j < std::size(M_); ++j){
      norm_n += std::pow(std::abs(M_[j][j]), n);
    }
    return std::pow(norm_n, 1./n);
  }


}//namespace
