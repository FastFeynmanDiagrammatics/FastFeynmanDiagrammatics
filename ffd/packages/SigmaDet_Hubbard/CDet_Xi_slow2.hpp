namespace ffd::sigmadet_hubbard{
  //  ffd::user_space::Timer t_temp, t_1, t_2, t_3, t_4;

  template<std::size_t n, class field_d = Real>
  void
  CDet_Xi_slow2_r(array2d<field_d, n, n> const& __restrict__ G0,
		  array3d<field_d, n, n, (1ul<<n)>&  __restrict__ Xi,
		  array1d<field_d, (1ul<<n)>& __restrict__ Z){
    static_assert(( n <= 16 ));
    //    t_3.ini();
    Xi.fill(0);
    
    Z[0] = 1;
    for( BinaryInt S = 1; S < (1ul<<n); ++S ){
      BinaryInt const compl_S = (1ul<<n)-1-S;
      int const card_S = ffd::set_theory::CardinalitySet(S);
      int const card_compl_S = n-card_S;
      bool const is_even_S = ((card_S&1) == 0);

      std::vector<BinaryInt> const digits_S = ffd::set_theory::VectorOfBinaryDigitsOf(S);
      vector2d<field_d> g0(card_S, card_S);
      //  t_2.ini();
      for( uint j = 0; j < card_S; ++j ){
	for( uint k = 0; k < card_S; ++k ){
	  g0(k, j) = G0(digits_S[k], digits_S[j]);
	}
      }
      // t_2.fin();
      //      t_temp.ini();
      auto const [I, det] = ffd::inverse_matrix_gauss::
        Inverse_and_Determinant(std::move(g0));
      //      t_temp.fin();
      auto const det2 = det*det;
      {
	if( is_even_S ){
	  Z[S] = det2;
	}else{
	  Z[S] = -det2;
	}
      }

      for( ulong j = 0; j < card_S; ++j ){
	for( ulong k = 0; k < card_S; ++k ){
	  auto const xi_temp = I[j*card_S+k]*det2;
	  Xi(digits_S[k], digits_S[j], S) = is_even_S ? -xi_temp : xi_temp;
	}
      }
    }

    // compute connected part
    //    t_1.ini();
    for(ulong V=1; V<(1ul<<n); ++V){
      auto const V_n = V*n*n;
      for(ulong S=((V-1)&V); S != 0; S=((S-1)&V) ){
	auto const S_n = S * n * n;
	auto const Z_V_S = Z[V-S];
	for( ulong j = 0; j < n; ++j ){
	  if( (S&(1ul<<j)) != 0 ){
	    auto const V_j_n = V_n + j * n;
	    auto const S_j_n = S_n + j * n;
	    std::array<field_d, n> temp;
	    temp.fill(Z_V_S);
#pragma unroll
	    for ( ulong k = 0; k < n; ++k ) {
	      temp[k] *= Xi[S_j_n+k];	      
	    } // for k in range(0, n)
	    // temp[0] *= Xi[S_j_n];
	    // if constexpr ( n > 1 ) {
	    // 	temp[1] *= Xi[S_j_n+1];
	    //   }
	    // if constexpr ( n > 2 ) {
	    // 	temp[2] *= Xi[S_j_n+2];
	    //   }
	    // if constexpr ( n > 3 ) {
	    // 	temp[3] *= Xi[S_j_n+3];
	    //   }
	    // if constexpr ( n > 4 ) {
	    // 	temp[4] *= Xi[S_j_n+4];
	    //   }
	    // if constexpr ( n > 5 ) {
	    // 	temp[5] *= Xi[S_j_n+5];
	    //   }
	    // if constexpr ( n > 6 ) {
	    // 	temp[6] *= Xi[S_j_n+6];
	    //   }
	    // if constexpr ( n > 7 ) {
	    // 	temp[7] *= Xi[S_j_n+7];
	    //   }
	    // if constexpr ( n > 8 ) {
	    // 	temp[8] *= Xi[S_j_n+8];
	    //   }
	    // if constexpr ( n > 9 ) {
	    // 	temp[9] *= Xi[S_j_n+9];
	    //   }
	    // if constexpr ( n > 10 ) {
	    // 	temp[10] *= Xi[S_j_n+10];
	    //   }
	    // if constexpr ( n > 11 ) {
	    // 	temp[11] *= Xi[S_j_n+11];
	    //   }
	    // if constexpr ( n > 12 ) {
	    // 	temp[12] *= Xi[S_j_n+12];
	    //   }
	    // if constexpr ( n > 13 ) {
	    // 	temp[13] *= Xi[S_j_n+13];
	    //   }
	    // if constexpr ( n > 14 ) {
	    // 	temp[14] *= Xi[S_j_n+14];
	    //   }
	    // if constexpr ( n > 15 ) {
	    // 	temp[15] *= Xi[S_j_n+15];
	    //   }
	    for ( ulong k = 0; k < n; ++k ) {
	      Xi[V_j_n+k] *= temp[k];
	    } // for j in range(0, n)
	  }
	}
      }
    }
    // t_1.fin();
    // t_3.fin();
    return;
  }

  
} // namespace
