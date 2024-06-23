namespace ffd::sigmadet_hubbard{

  template<typename nilpoly_t = ffd::nilpotent_polynomial::NilpotentPolynomial<Real>,
	   typename G0_t = int>
  
  auto
  compute_Xi(G0_t const& G0){
    using full_nilpoly_t = decltype(G0[0].extend_to_set(0, 0));
    int const n = ffd::core_math::sqrt_int( size(G0) );
    BinaryInt const two_n = (1<<n);
    full_nilpoly_t Z(n);
    
    
    Z[0] = 1.;
    std::vector<full_nilpoly_t> Xi(size(G0), full_nilpoly_t(n));
    // for(uint j=0; j<n*n; ++j){
    //   Xi[j] = nilpoly_t(n);
    // }
    for( BinaryInt S = 1; S < two_n; ++S ){
      BinaryInt const compl_S = two_n-1-S;
      int const card_S = ffd::set_theory::CardinalitySet(S);
      int const card_compl_S = n-card_S;
      bool const is_even_S = ((card_S&1) == 0);
      std::vector<BinaryInt> const digits_S = ffd::set_theory::VectorOfBinaryDigitsOf(S);
      std::vector<nilpoly_t> g0(card_S*card_S, card_compl_S);
      for( int j = 0; j < card_S; ++j ){
	for( int k = 0; k < card_S; ++k ){
	  g0[j*card_S + k] = G0[digits_S[j]*n + digits_S[k]].restrict_to_set(compl_S);
	}
      }


      auto const [I, det] = ffd::inverse_matrix_gauss::
	Inverse_and_Determinant(std::move(g0));
      auto const det2 = det*det;
      {
	auto const det2_Z = det2.extend_to_set(compl_S, n).shift(S, n);
	// auto const det2_Z = det2.extend_to_set(compl_S).shift(S);
	if( is_even_S ){
	  Z += det2_Z;
	}else{
	  Z -= det2_Z;
	}
      }
      
      
      for( int j = 0; j < card_S; ++j ){
	for( int k = 0; k < card_S; ++k ){
	  auto const xi_temp = (I[j*card_S+k]*det2).extend_to_set(compl_S, n).shift(S, n);
	  // auto const xi_temp = (I[j*card_S+k]*det2).extend_to_set(compl_S).shift(S);
	  if( !is_even_S ){
	    Xi[digits_S[j]*n+digits_S[k]] += xi_temp;
	  }else{
	    Xi[digits_S[j]*n+digits_S[k]] -= xi_temp;
	  }
	}
      }
    }
    
    
    for( int j = 0; j < n; ++j ){
      for( int k = 0; k < n; ++k ){
    	  Xi[j*n+k] /= Z;
      }
    }
    
    return Xi;
  }

}//namespace
