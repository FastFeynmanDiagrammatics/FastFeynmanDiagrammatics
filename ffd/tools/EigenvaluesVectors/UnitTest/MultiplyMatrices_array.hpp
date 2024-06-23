

namespace ffd::eigenvalues_vectors::unit_test{
  
  using namespace std;


  template<typename Field, auto n>
  auto MultiplyMatrices_array(std::array<std::array<Field, n>, n> const& M1_,
			      std::array<std::array<Field, n>, n> const& M2_){
    std::array<std::array<Field, n>, n> ret;
    for(decltype(n) j=0; j<n; ++j){
      for(decltype(n) k=0; k<n; ++k){
	ret[j][k] = 0.;
	for(decltype(n) m=0; m<n; ++m){
	  ret[j][k] += M1_[j][m]*M2_[m][k];
	}
      }
    }
    return ret;
  }


}//namespace
