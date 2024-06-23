namespace ffd::core_math{

  template<int n, int k>
  
  constexpr int
  BinomialCoefficient(){
    static_assert(n >= 0);
    static_assert(k >= 0);
    if constexpr(k>n){
	return 0;
      }
    if constexpr(n==0 || k==0 || k==n){
	return 1;
      }else{
      return BinomialCoefficient<n-1, k-1>()+BinomialCoefficient<n-1, k>();
    }
  }


  constexpr int
  BinomialCoefficient(int n, int k){
    assert(n >= 0);
    assert(k >= 0);
    if(k>n){
	return 0;
      }
    if(n==0 || k==0 || k==n){
	return 1;
      }else{
      return BinomialCoefficient(n-1, k-1)+BinomialCoefficient(n-1, k);
    }
  }


}//namespace
