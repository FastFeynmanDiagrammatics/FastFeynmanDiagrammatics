namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  [[nodiscard]] NilpotentPolynomial<Field, SetType>
  bound_divide(const NilpotentPolynomial<Field, SetType>& P,
               const NilpotentPolynomial<Field, SetType>& Q,
               const BinaryInt max_order) noexcept{
    assert(( Q.NotZero ));

    if( !P.NotZero ){
      return NilpotentPolynomial<Field, SetType>();
    }
    else if(Q.size() == 1){
      return P/Q[0];
    }

    NilpotentPolynomial<Field, SetType> ret( int(Q.linear_size()) );
    std::size_t const P_size = P.size();

    Field const one_over_Q0 = 1./Q[0];
    for(std::size_t set = 0ul; set < P_size; ++set){
      ret[set] = P[set] * one_over_Q0;
    }

    NilpotentPolynomial<Field, SetType> q = Q*one_over_Q0;

    Field const ret0 = ret[0];

    if constexpr (SetType == Full){
      for(BinaryInt V = 1u; V < P_size; ++V){
        if (__builtin_popcount(V)<=max_order){
          ret[V] -= q[V] * ret0;
          BinaryInt const V_max = V-(1<<(31-__builtin_clz(V)));
          for(BinaryInt S = ((V-1)&V); S > V_max; S = ((S-1)&V)){
            BinaryInt const V_S = V-S;
            ret[V] -= q[S] * ret[V_S] + q[V_S] * ret[S];
          }
        }
      }
    }
    else if constexpr (SetType == Even){
      for(BinaryInt V = 1u; V < P_size; ++V){
        if (!__builtin_parity(V) && __builtin_popcount(V)<=max_order){
          ret[V] -= q[V] * ret0;
          BinaryInt const V_max = V-(1<<(31-__builtin_clz(V)));
          for(BinaryInt S = ((V-1)&V); S > V_max; S = ((S-1)&V)){
            if (!__builtin_parity(S))
            {
              BinaryInt const V_S = V-S;
              ret[V] -= q[S] * ret[V_S] + q[V_S] * ret[S];
            }
          }
        }
      }
    }
    return ret;
  }

}//namespace
