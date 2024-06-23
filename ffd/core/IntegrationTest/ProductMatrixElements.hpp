

namespace ffd::core_integration_test{

  template<typename Field, int intType>
  Field ProductMatrixElements(ffd::wick_matrix::WickMatrix<Field, intType> const& M_,
			      std::vector<char> const& permutation_){
    Field ret = 1.;
    using std::size;
    for(char j=0; j < size(M_); ++j){
      ret *= M_(j, permutation[j]);
    }
    return ret;
  }
  
  
}//namespace 
