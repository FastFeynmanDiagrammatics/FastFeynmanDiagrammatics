namespace ffd::laurent_polynomial {

  template < class F, int m, int M >
  bool
  operator< (
	     P<F, m, M> const& P0,
	     P<F, m, M> const& P1
	     ) {
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      if ( P0.coef[j] < P1.coef[j] ) return true;
      else if ( P0.coef[j] > P1.coef[j] ) return false;
    } // for j in range(0, -m+M+1)
    return false;
  }

  template < class F, int m, int M >
  bool
  operator== (
	     P<F, m, M> const& P0,
	     P<F, m, M> const& P1
	     ) {
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      if ( P0.coef[j] != P1.coef[j] ) return false;
    } // for j in range(0, -m+M+1)
    return true;
  }

  template < class F, int m, int M >
  bool
  operator!= (
	     P<F, m, M> const& P0,
	     P<F, m, M> const& P1
	     ) {
    return !(P0 == P1);
  }

  template < class F, int m, int M >
  bool
  operator> (
	     P<F, m, M> const& P0,
	     P<F, m, M> const& P1
	     ) {
    return P1 < P0;
  }

  template < class F, int m, int M >
  bool
  operator< (
	     P<F, m, M> const& P0,
	     F x
	     ) {
    Real sum = 0;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      sum += P0.coef[j];
    } // for j in range(0, -m+M+1)
    return sum < x;
  }
  
} // namespace ffd::laurent_polynomial
