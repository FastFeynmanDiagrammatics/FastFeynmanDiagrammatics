namespace ffd::sigmadet_hubbard{

  template<int n,
	   typename G0_t = int>
  
  auto
  compute_Xi_CDet_slow(G0_t const& G0){
    using nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real>;
    auto const two_n = (1u<<n);
    nilpoly_t Z(n);
    
    Z[0] = 1.;
    array2d<nilpoly_t, n, n> Xi;
    Xi.fill(nilpoly_t(n));
    for ( ulong j = 0; j < n; ++j ) {
      for ( ulong k = 0; k < n; ++k ) {
	Xi(j, k)[0] = 0.;
      } // for k in range(0, n)
    } // for j in range(0, n)
    for( BinaryInt S = 1; S < two_n; ++S ){
      BinaryInt const compl_S = two_n-1-S;
      int  const card_S = ffd::set_theory::CardinalitySet(S);
      bool const is_even_S = ((card_S&1) == 0);
      std::vector<BinaryInt> const digits_S = ffd::set_theory::VectorOfBinaryDigitsOf(S);
      vector2d<Real> g0 (card_S, card_S, 0.);
      for( int j = 0; j < card_S; ++j ){
	for( int k = 0; k < card_S; ++k ){
	  g0(k, j) = G0[digits_S[j]*n + digits_S[k]];
	}
      }
      

      auto const [I, det] = ffd::inverse_matrix_gauss::
	Inverse_and_Determinant(std::move(g0));
      // std::cerr << S << " " << det;
      // for ( ulong j = 0; j < n; ++j ) {
      // 	for ( ulong k = 0; k < n; ++k ) {
      // 	  std::cerr << I[j+k*n] << " ";
      // 	} // for k in range(0, n)
      // } // for j in range(0, n)
      // std::cerr << "\n\n";
      auto const det2 = det*det;
      {
	if( is_even_S ){
	  Z[S] += det2;
	}else{
	  Z[S] -= det2;
	}
      }
      
      
      for( int j = 0; j < card_S; ++j ){
	for( int k = 0; k < card_S; ++k ){
	  auto const xi_temp = I[j*card_S+k]*det2;
	  if( !is_even_S ){
	    Xi[digits_S[j]*n+digits_S[k]][S] += xi_temp;
	  }else{
	    Xi[digits_S[j]*n+digits_S[k]][S] -= xi_temp;
	  }
	}
      }
    }
    
    
    for( int j = 0; j < n; ++j ){
      for( int k = 0; k < n; ++k ){
	Xi[j*n+k] /= Z;
      }
    }
    
    return std::make_pair(Xi, Z);
  }

  
} // namespace
