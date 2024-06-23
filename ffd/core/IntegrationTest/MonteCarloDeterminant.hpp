

namespace ffd::core_integration_test{

  template<typename Field, typename IntType>
  Field MonteCarloDeterminant(ffd::wick_matrix::WickMatrix<Field, IntType> const& M_,
			      long NumSamples_, long NumShuffles_=0){
    using std::size;
    std::vector<char> permutation(size(M_));
    for(char j=0; j < size(M_); ++j){
      permutation[j] = j;
    }
    Field ret = 0.l;
    int sign = 1;
    for(long t=0; t < NumSamples_; ++t, sign *= -1){
      for(long u=0; u < 2*NumShuffles_ + 1; ++u){
	ShufflePermutation(permutation);
      }
      Field factor = sign;
      for(int j=0; j < size(M_); ++j){
	factor *= M_(j, permutation[j]);
      }
      ret += factor;
    }
    return ffd::core_math::Factorial(size(M_))*ret / NumSamples_;
  }
  
}//namespace
