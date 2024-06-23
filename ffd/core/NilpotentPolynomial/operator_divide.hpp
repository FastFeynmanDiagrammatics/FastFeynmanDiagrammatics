namespace ffd::nilpotent_polynomial {

template <typename Field, SetT SetType>
[[nodiscard]] NilpotentPolynomial<Field, SetType>
NilpotentPolynomial<Field, SetType>::operator/(
    const NilpotentPolynomial<Field, SetType>& Q) const noexcept {
  assert((Q.NotZero));

  // ADD [[unlikely]]
  if (!this->NotZero) {
    return NilpotentPolynomial<Field, SetType>();
  } else if (Q.size() == 1) {
    return *this / Q[0];
  }

  NilpotentPolynomial<Field, SetType> ret(int(Q.linear_size()));
  std::size_t const this_size = this->size();
  Field const one_over_Q0 = 1. / Q[0];
  for (std::size_t set = 0ul; set < this_size; ++set) {
    ret[set] = (*this)[set] * one_over_Q0;
  }

  NilpotentPolynomial<Field, SetType> q = Q * one_over_Q0;

  Field const ret0 = ret[0];

  if constexpr (SetType == Full) {
    for (BinaryInt V = 1u; V < this_size; ++V) {
      ret[V] -= q[V] * ret0;
      BinaryInt const V_max = V - (1 << (31 - __builtin_clz(V)));
      for (BinaryInt S = ((V - 1) & V); S > V_max; S = ((S - 1) & V)) {
        BinaryInt const V_S = V - S;
        ret[V] -= q[S] * ret[V_S] + q[V_S] * ret[S];
      }
    }

    // // OLD
    // for(BinaryInt V = 1u; V < this_size; ++V){
    //   for(BinaryInt S = V; S != 0u; S = ((S-1)&V)){
    // 		ret[V] -= q[S] * ret[V-S];
    // 	}
    // }

  } else if constexpr (SetType == Even) {
    for (BinaryInt V = 1u; V < this_size; ++V) {
      if (!__builtin_parity(V)) {
        ret[V] -= q[V] * ret0;
        BinaryInt const V_max = V - (1 << (31 - __builtin_clz(V)));
        for (BinaryInt S = ((V - 1) & V); S > V_max; S = ((S - 1) & V)) {
          if (!__builtin_parity(S)) {
            BinaryInt const V_S = V - S;
            ret[V] -= q[S] * ret[V_S] + q[V_S] * ret[S];
          }
        }
      }
    }
  }
  // else if constexpr (SetType == Odd){
  //   for(BinaryInt V = 1u; V < this_size; ++V){
  //     if (__builtin_parity(V)){
  //       ret[V] -= q[V] * ret[0];
  //       BinaryInt V_max = V-(1<<(31-__builtin_clz(V)));
  //       for(BinaryInt S = ((V-1)&V); S > V_max; S = ((S-1)&V)){
  //         if (__builtin_parity(S)){
  //           BinaryInt V_S = V-S;
  //           ret[V] -= q[S] * ret[V_S];
  //           ret[V] -= q[V_S] * ret[S];
  //         }
  //       }
  //     }
  //   }
  // }
  return ret;
}

template <typename Field, SetT SetType>
NilpotentPolynomial<Field, SetType>&
NilpotentPolynomial<Field, SetType>::operator/=(
    const NilpotentPolynomial<Field, SetType>& den_) noexcept {
  return *this = *this / den_;
}

}  // namespace ffd::nilpotent_polynomial
