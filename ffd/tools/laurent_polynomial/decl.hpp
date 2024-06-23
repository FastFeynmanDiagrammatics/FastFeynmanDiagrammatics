namespace ffd::laurent_polynomial {
  
  using ulong = unsigned long;
  
  template < class Field,
	     int m, int M >
  struct P {
    static_assert ( m <= 0 );
    static_assert ( M >= 0 );
    
    std::array<Field, -m+M+1> coef;
    
    P() {}
    template < class F,
	       std::enable_if_t<std::is_same_v<Field, F> ||
				std::is_same_v<F, int>, int> = 0 > 
    P(F const& x) {coef.fill(Field(0)); coef[-m] = x;}

    inline Field const& operator[](int j) const { return coef[-m+j];}
    inline Field& operator[](int j){ return coef[-m+j];}
  };

  template < class F, int m, int M >
  P<F, m, M>
  operator+(
	    P<F, m, M> const& P0
	    ) {
    return P0;
  }

  template < class F, int m, int M >
  P<F, m, M>
  operator-(
	    P<F, m, M> const& P0
	    ) {
    P<F, m, M> ret;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret.coef[j] = - P0.coef[j];
    } // for j in range(0, -m+M+1)
    return ret;
  }
  
  template < class F, int m, int M >
  P<F, m, M>
  operator+(
	    P<F, m, M> const& P0,
	    P<F, m, M> const& P1
	    ) {
    P<F, m, M> ret = P0;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret.coef[j] += P1.coef[j];
    } // for j in range(0, -m+M+1)
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>
  operator+(
	    P<F, m, M> const& P0,
	    F x
	    ) {
    auto ret = P0;
    ret.coef[-m] += x;
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>
  operator+(
	    F x,
	    P<F, m, M> const& P0
	    ) {
    auto ret = P0;
    ret.coef[-m] += x;
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>&
  operator+=(
	     P<F, m, M>      & P0,
	     P<F, m, M> const& P1
	     ) {
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      P0.coef[j] += P1.coef[j];
    } // for j in range(0, -m+M+1)
    return P0;
  }

  template < class F, int m, int M >
  P<F, m, M>&
  operator+=(
	     P<F, m, M>      & P0,
	     F x
	     ) {
    P0.coef[-m] += x;
    return P0;
  }
  
  template < class F, int m, int M >
  P<F, m, M>
  operator-(
	    P<F, m, M> const& P0,
	    P<F, m, M> const& P1
	    ) {
    P<F, m, M> ret = P0;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret.coef[j] -= P1.coef[j];
    } // for j in range(0, -m+M+1)
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>
  operator-(
	    P<F, m, M> const& P0,
	    F x
	    ) {
    auto ret = P0;
    ret.coef[-m] -= x;
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>
  operator-(
	    F x,
	    P<F, m, M> const& P0
	    ) {
    auto ret = -P0;
    ret.coef[-m] += x;
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>&
  operator-=(
	     P<F, m, M>      & P0,
	     P<F, m, M> const& P1
	     ) {
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      P0.coef[j] -= P1.coef[j];
    } // for j in range(0, -m+M+1)
    return P0;
  }

  
  template < class F, int m, int M >
  P<F, m, M>&
  operator-=(
	     P<F, m, M>      & P0,
	     F x
	     ) {
    P0.coef[-m] -= x;
    return P0;
  }

  
  template < class F, int m, int M >
  P<F, m, M>
  operator*(
	    P<F, m, M> const& P0,
	    P<F, m, M> const& P1
	    ) {
    P<F, m, M> ret;
    if constexpr ( m == -1 && M == 1 ) {
	ret.coef[0] = P0.coef[0] * P1.coef[1] + P0.coef[1] * P1.coef[0];
	ret.coef[1] = P0.coef[0] * P1.coef[2] + P0.coef[1] * P1.coef[1] + P0.coef[2] * P1.coef[0];
	ret.coef[2] = P0.coef[1] * P1.coef[2] + P0.coef[2] * P1.coef[1];
      }
    if constexpr ( m == -2 && M == 2) {
	ret[-2] = P0[-2]*P1[0]+P0[-1]*P1[-1]+P0[0]*P1[-2];
	ret[-1] = P0[-2]*P1[1]+P0[-1]*P1[0]+P0[0]*P1[-1]+P0[1]*P1[-2];
	ret[0] = P0[-2]*P1[2]+P0[-1]*P1[1]+P0[0]*P1[0]+P0[1]*P1[-1]+P0[2]*P1[-2];
	ret[1] = P0[-1]*P1[2]+P0[0]*P1[1]+P0[1]*P1[0]+P0[2]*P1[-1];
	ret[2] = P0[0]*P1[2]+P0[1]*P1[1]+P0[2]*P1[0];
      }
    return ret;
  }

  
  template < class F, int m, int M >
  P<F, m, M>
  operator*(
	    P<F, m, M> const& P0,
	    F x
	    ) {
    P<F, m, M> ret;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret.coef[j] = P0.coef[j]*x;
    } // for j in range(0, -m+M+1)
    return ret;
  }

  
  template < class F, int m, int M >
  P<F, m, M>
  operator*(
	    F x,
	    P<F, m, M> const& P0
	    ) {
    P<F, m, M> ret;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret.coef[j] = P0.coef[j]*x;
    } // for j in range(0, -m+M+1)
    return ret;
  }

  
  template < class F, int m, int M >
  P<F, m, M>&
  operator*=(
	     P<F, m, M>      & P0,
	     P<F, m, M> const& P1
	    ) {
    return P0 = P0 * P1;
  }

  
  template < class F, int m, int M >
  P<F, m, M>&
  operator*=(
	     P<F, m, M>      & P0,
	     F x
	    ) {
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      P0.coef[j] *= x;
    } // for j in range(0, -m+M+1)
    return P0;
  }

  
  template < class F, int m, int M >
  P<F, m, M>
  operator/(
	    P<F, m, M> const& P0,
	    F x
	    ) {
    return P0 * (1/x);
  }
  
  template < class F0, int m, int M >
  auto
  abs_norm ( P<F0, m, M> const& P0 ) {
    using std::abs;
    std::decay_t<decltype(abs(P0[0]))> ret = 0;
    for ( ulong j = 0; j < -m+M+1; ++j ) {
      ret += abs(P0.coef[j]);
    } // for j in range(0, -m+M+1)
    return ret;
  }

  template < class F, int m, int M >
  int
  leading_order (
		 P<F, m, M> const& P0
		 ) {
    auto const eps = 1e-12;
    using std::abs;
    auto norm_abs_eps = abs_norm(P0)*eps;
    for ( int j = 0; j < -m+M+1; ++j ) {
      if ( abs(P0.coef[j]) > norm_abs_eps ) {
	return m+j; 
      }
    } // for j in range(0, -m+M+1)
    return M;
  }
  
  template < class Fx, class F0, int m, int M,
	     std::enable_if_t<std::is_same_v<Fx, F0> || std::is_same_v<Fx, int>, int> = 0>
  P<F0, m, M>
  operator/(
	    Fx x,
	    P<F0, m, M> const& P0
	    ) {
    P<F0, m, M> ret;
    if constexpr ( m == -1 && M == 1) {
	switch(leading_order(P0)){
	case -1:
	  ret.coef[0] = 0;
	  ret.coef[1] = 0;
	  ret.coef[2] = x / P0.coef[0];
	  break;
	case 0:
	  ret.coef[0] = 0;
	  ret.coef[1] = x / P0.coef[1];
	  ret.coef[2] = - x * P0.coef[2] / (P0.coef[1]*P0.coef[1]);
	  break;
 	case 1:
	  ret.coef[0] = x / P0.coef[2];
	  ret.coef[1] = 0;
	  ret.coef[2] = 0;
	  break;
	}
      }
    if constexpr ( m == -2 && M == 2 ) {
	switch(leading_order(P0)){
	case -2:
	  ret[-2] = 0;
	  ret[-1] = 0;
	  ret[0] = 0;
	  ret[1] = 0;
	  ret[2] = x / P0[-2];
	  break;
	case -1:
	  ret[-2] = 0;
	  ret[-1] = 0;
	  ret[0] = 0;
	  ret[1] = x / P0[-1];
	  ret[2] = -x*P0[0] / (P0[-1]*P0[-1]);
	  break;
	case 0:
	  ret[-2] = 0;
	  ret[-1] = 0;
	  ret[0] = x / P0[0];
	  ret[1] = -x * P0[1] / (P0[0]*P0[0]);
	  ret[2] = -x * P0[2] / (P0[0]*P0[0])+ x*P0[1]*P0[1]/(P0[0]*P0[0]*P0[0]);
	case 1:
	  ret[-2] = 0;
	  ret[-1] = x/P0[1];
	  ret[0] = -x*P0[2]/(P0[1]*P0[1]);
	  ret[1] = x*P0[2]*P0[2] / (P0[1]*P0[1]*P0[1]);
	  ret[2] = -x*P0[2]*P0[2]*P0[2] / (P0[1]*P0[1]*P0[1]);
	  break;
	case 2:
	  ret[-2] = x / P0[2];
	  ret[-1] = 0;
	  ret[0] = 0;
	  ret[1] = 0;
	  ret[2] = 0;
	  break;
	}
      }
    
    return ret;
  }  
  
  template < class F, int m, int M >
  P<F, m, M>
  operator/(
	    P<F, m, M> const& P0,
	    P<F, m, M> const& P1
	    ) {
    P<F, m, M> ret;
    if constexpr ( m == -1 && M == 1) {
	switch(leading_order(P1)){
	case -1:
	  ret.coef[0] = 0;
	  ret.coef[1] = P0.coef[0] / P1.coef[0];
	  ret.coef[2] = P0.coef[1] / P1.coef[0] - P0.coef[0]*P1.coef[1] / (P1.coef[0]*P1.coef[0]);
	  break;
	case 0:
	  ret.coef[0] = P0.coef[0] / P1.coef[1];
	  ret.coef[1] = P0.coef[1] / P1.coef[1] - P0.coef[0]*P1.coef[2] / (P1.coef[1]*P1.coef[1]);
	  ret.coef[2] = P0.coef[2] / P1.coef[1]
	    -P0.coef[1]*P1.coef[2]/(P1.coef[1]*P1.coef[1])
	    +P0.coef[0]*P1.coef[2]*P1.coef[2] / (P1.coef[1]*P1.coef[1]*P1.coef[1]);
	  break;
	case 1:
	  ret.coef[0] = P0.coef[1] / P1.coef[2];
	  ret.coef[1] = P0.coef[2] / P1.coef[2];
	  ret.coef[2] = 0;
	  break;
	}
      }
    if constexpr ( m == -2 && M == 2) {
	ret = P0*(1/P1);
      }
    return ret;
  }

  template < class F, int m, int M >
  P<F, m, M>&
  operator/=(
	     P<F, m, M>      & P0,
	     P<F, m, M> const& P1
	     ) {
    return P0 = P0 * P1;
  }
  
  template < class F, int m, int M >
  P<F, m, M>&
  operator/=(
	     P<F, m, M> & P0,
	     F x
	     ) {
    return P0 *= 1/x;
  }
  
} // namespace ffd::laurent_polynomial
