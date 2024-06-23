namespace ffd::laurent_polynomial {

  template < class F, int m, int M >
  auto
  to_sstream (P<F, m, M> const& P0) {
    static_assert(m==-1);
    static_assert(M==1);
    std::stringstream ss;
    if ( P)
  }
  
} // namespace ffd::laurent_polynomial
