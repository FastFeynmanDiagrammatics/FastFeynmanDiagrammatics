

namespace ffd::conformal_mappings{

  template<typename Field>
  Real ZinnJustin2<Field>::CoefficientMapping(int n, int m){
    assert(n >= 1);
    assert(n <= m);
    long double ret = pow(4., n);
    for(int k=0; k < m-n; ++k){
      ret *= 2*n+k;
      if(k >= 2){
	ret /= k;
      }
    }
    return ret;
  }
  

}//namespace
