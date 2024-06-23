namespace ffd::four_point_vertex_hubbard{

  template<int order,
	   typename G0_t,
	   typename nilpoly_t>
  
  auto
  connected(G0_t const& G0){
    int const n = ffd::core_math::sqrt_int( size(G0) );
    assert((n==order));
    BinaryInt constexpr two_order = (1<<order);
    nilpoly_t Z(order);
    std::array<nilpoly_t, order*order*order*order> chi;
    chi.fill(nilpoly_t(order));
    
    
    Z[0] = 1.;
    for( BinaryInt S = 1; S < two_order; ++S ){
      BinaryInt const compl_S = two_order-1-S;
      int const card_S = ffd::set_theory::
	CardinalitySet(S);
      int const card_compl_S = order - card_S;
      bool const is_even_S = ((card_S&1) == 0);
      std::vector<BinaryInt> const digits_S =
	ffd::set_theory::VectorOfBinaryDigitsOf(S);
      std::vector<nilpoly_t> g0(card_S*card_S, card_compl_S);
      for( int j = 0; j < card_S; ++j ){
	for( int k = 0; k < card_S; ++k ){
	  g0[j*card_S + k] = G0[digits_S[j]*order + digits_S[k]].restrict_to_set(compl_S);
	}
      }
      
      
      auto const [I, det] = ffd::inverse_matrix_gauss::
	Inverse_and_Determinant(std::move(g0));
      auto const det2 = det*det;
      {
	auto const det2_Z = det2.extend_to_set(compl_S, order).shift(S, order);
	// auto const det2_Z = det2.extend_to_set(compl_S).shift(S);
	if( is_even_S ){
	  Z += det2_Z;
	}else{
	  Z -= det2_Z;
	}
      }
      
      
      for( int j1 = 0; j1 < card_S; ++j1 ){
	for( int k1 = 0; k1 < card_S; ++k1 ){
	  for( int j2 = 0; j2 < card_S; ++j2 ){
	    for( int k2 = 0; k2 < card_S; ++k2 ){
	      auto const chi_temp = (I[j1*card_S+k1]*
				     I[j2*card_S+k2]*
				     det2).
		extend_to_set(compl_S, n).shift(S, n);
	      auto& chi_r = chi[((digits_S[j1]*order+digits_S[k1])*order+
				 digits_S[j2])*order+digits_S[k2]];
	      if( !is_even_S ){
		chi_r += chi_temp;
	      }else{
		chi_r -= chi_temp;
	      }
	    }
	  }
	}
      }
    }
    
    
    for( int j1 = 0; j1 < order; ++j1 ){
      for( int k1 = 0; k1 < order; ++k1 ){
	for( int j2 = 0; j2 < order; ++j2 ){
	  for( int k2 = 0; k2 < order; ++k2 ){
	    chi[((j1*order+k1)*order+j2)*order+k2] /= Z;
	  }
	}
      }
    }
    
    
    return chi;
  }
    
}//namespace
