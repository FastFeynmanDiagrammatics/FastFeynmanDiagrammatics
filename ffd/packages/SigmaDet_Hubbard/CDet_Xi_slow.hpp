namespace ffd::sigmadet_hubbard{
  //  ffd::user_space::Timer t_temp, t_1, t_2, t_3, t_4;

  template<std::size_t n, class field_d = Real>
  void
  CDet_Xi_slow_r(array2d<field_d, n, n> const& __restrict__ G0,
		 array3d<field_d, n, n, (1ul<<n)>&  __restrict__ Xi) {
    //    t_3.ini();
    Xi.fill(0);
    array1d<field_d, (1ul<<n)> Z;
    
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
	auto const V_S = V - S;
	for( ulong j = 0; j < n; ++j ){
	  if( (S&(1ul<<j)) != 0 ){
	    auto const V_j_n = V_n + j * n;
	    auto const S_j_n = S_n + j * n;
#pragma unroll
	    for( ulong k = 0; k < n; ++k ){
	      Xi[V_j_n+k] -= Xi[S_j_n+k]*Z[V_S];
	    }
	  }
	}
      }
    }
    // t_1.fin();
    // t_3.fin();
    return;
  }


  template<std::size_t n,
  	   class field_d = Real>

  void
  CDet_Xi_slow_r(array2d<field_d, n, n> const& G0,
		 array3d<field_d, (1ul<<n), n, n>& Xi,
		 array1d<field_d, (1ul<<n)>& Z){
    Xi.fill(0);
    
    Z[0] = 1.;
    for( BinaryInt S = 1; S < (1ul<<n); ++S ){
      BinaryInt const compl_S = (1ul<<n)-1-S;
      int const card_S = ffd::set_theory::CardinalitySet(S);
      int const card_compl_S = n-card_S;
      bool const is_even_S = ((card_S&1) == 0);
      std::vector<BinaryInt> const digits_S = ffd::set_theory::VectorOfBinaryDigitsOf(S);
      vector2d<field_d> g0(card_S, card_S);
      for( int j = 0; j < card_S; ++j ){
  	for( int k = 0; k < card_S; ++k ){
  	  g0(k, j) = G0(digits_S[k], digits_S[j]);
  	}
      }
      
      
      auto const [I, det] = ffd::inverse_matrix_gauss::
  	Inverse_and_Determinant(std::move(g0));
      auto const det2 = det*det;
      {
  	if( is_even_S ){
  	  Z[S] = det2;
  	}else{
  	  Z[S] = -det2;
  	}
      }
      
      
      for( int j = 0; j < card_S; ++j ){
  	for( int k = 0; k < card_S; ++k ){
  	  auto const xi_temp = I[j*card_S+k]*det2;
  	  if( !is_even_S ){
  	    Xi(S, digits_S[k], digits_S[j]) = xi_temp;
  	  }else{
	    Xi(S, digits_S[k], digits_S[j]) = -xi_temp;
  	  }
  	}
      }
    }

    // compute connected part
    for( int j = 0; j < n; ++j ){
      for( int k = 0; k < n; ++k ){
  	for(uint V=1; V<(1ul<<n); ++V){
  	  for(uint S=((V-1)&V); S>0; S=((S-1)&V) ){
  	    Xi(V, k, j) -= Xi(S, k, j) * Z[V-S];
  	  }
  	}
      }
    }

    
    return;
  }

}//namespace
