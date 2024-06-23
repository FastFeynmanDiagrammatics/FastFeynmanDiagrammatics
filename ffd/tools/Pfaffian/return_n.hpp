namespace ffd::gauss_pfaffian{

  inline int size_n(int n){ return n*(2*n-1); }

  
  template<typename size_t>
  int return_n(size_t size_v){
    assert(( size_v >= 0 ));
    int n = 0;
    for(; size_v != size_n(n); ++n){}
    return n;
  }


}//namespace
