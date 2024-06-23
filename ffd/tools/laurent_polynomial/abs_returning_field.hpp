//namespace ffd::laurent_polynomial::abs_returning_field {
namespace ffd::laurent_polynomial {

  template < class F, int m, int M >
  P<F, m, M>
  abs (P<F, m, M> const& P0) {
    using std::abs;
    P<F, m, M> ret;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret.coef[j] = abs(P0.coef[j]); 
    } // for j in range(0, -m+M+1)
    return ret;
  }
  
} // namespace ffd::laurent_polynomial::abs_returning_field
