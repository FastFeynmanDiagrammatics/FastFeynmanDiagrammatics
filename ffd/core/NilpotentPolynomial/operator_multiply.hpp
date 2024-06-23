namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>

  NilpotentPolynomial<Field, SetType>

  NilpotentPolynomial<Field, SetType>::
  operator*(NilpotentPolynomial<Field, SetType> const& P) const noexcept{
    if(!this->NotZero || !P.NotZero){
      return NilpotentPolynomial<Field, SetType>();
    }


    std::size_t this_size = this->size();
    if(this_size == 1){
      return P * ((*this)[0]);
    }else if(P.size() == 1){
      return *this * P[0];
    }

    NilpotentPolynomial<Field, SetType> ret( int(this->linear_size()) );

    if constexpr (SetType == Full){
      ret[0] += P[0] * this->coef[0];
      for(BinaryInt V = 1u; V < this_size; ++V){
        BinaryInt const V_max = V-(1<<(31-__builtin_clz(V)));
        for(BinaryInt S = V; S > V_max; S = ((S-1)&V)){
          BinaryInt const V_S = V-S;
          ret[V] += P[S] * this->coef[V_S] + P[V_S] * this->coef[S];
        }
      }

		// // OLD
		// NilpotentPolynomial<Field, SetType> ret = *this*P[0];
		// for(BinaryInt V = 1u; V < this_size; ++V){
		//   for(BinaryInt S = V; S != 0; S = ((S-1)&V)){
		// 	  ret[V] += P[S] * this->coef[V-S];
		//   }
		// }
	}
	else if constexpr (SetType == Even){
	    ret[0] += P[0] * this->coef[0];
	    for(BinaryInt V = 1u; V < this_size; ++V){
	      if (!__builtin_parity(V)){
	        BinaryInt const V_max = V-(1<<(31-__builtin_clz(V)));
	        for(BinaryInt S = V; S > V_max; S = ((S-1)&V)){
	          if (!__builtin_parity(S))
		    {
	            BinaryInt const V_S = V-S;
	            ret[V] += P[S] * this->coef[V_S] + P[V_S] * this->coef[S];
	          }
	        }
	      }
	    }
	}
  // else if constexpr (SetType == Odd){
  //   ret[0] += P[0] * this->coef[0];
  //   for(BinaryInt V = 1u; V < this_size; ++V){
  //     if (__builtin_parity(V)){
  //       BinaryInt V_max = V-(1<<(31-__builtin_clz(V)));
  //       for(BinaryInt S = V; S > V_max; S = ((S-1)&V)){
  //         if (__builtin_parity(S)){
  //           BinaryInt V_S = V-S;
  //           ret[V] += P[S] * this->coef[V_S] + P[V_S] * this->coef[S];
  //         }
  //       }
  //     }
  //   }
  // }

    return ret;
  }



  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>&
  NilpotentPolynomial<Field, SetType>::operator*=(const NilpotentPolynomial<Field, SetType>& P) noexcept{
    return *this = *this * P;
  }


}//namespace ffd::nilpotent_polynomial
