namespace ffd::heat_bath_mc{

  template<int max_order,
	   bool is_translation_invariant,
	   typename seed_t,
	   typename vector_coord_t>
  
  auto
  create_new_vertices(seed_t const& Phi,
		      vector_coord_t const& X,
		      BinaryInt S){
    using element_t = typename ffd::type_traits::element_of_container_t<vector_coord_t>;
    std::array<element_t, max_order> X_new;
    std::vector<BinaryInt> const digits_S =
      ffd::set_theory::VectorOfBinaryDigitsOf(S);
    BinaryInt const card_S = __builtin_popcount(S);
    int const num_vertices_to_add = max_order-int(card_S);
    
    
    for(BinaryInt j=0; j<card_S; ++j){
      X_new[j] = X[digits_S[j]];
    }
    
    
    for(int j=0; j<num_vertices_to_add; ++j){
      if ( !is_translation_invariant || !(card_S == 0 && j == 0) ) {
	int vertex = ffd::random_distributions::
	  RandomInRange(-(!is_translation_invariant), card_S+j);
	if constexpr (!is_translation_invariant){
	    if (vertex<0) {
	      X_new[card_S+j] = Phi(X[max_order]);	    
	    } else {
	      X_new[card_S+j] = Phi(X_new[vertex]);	    
	    }
	  } else {
	  X_new[card_S+j] = Phi(X_new[vertex]);
	}
      } else {
	X_new[0] = X[0];
      }
    }
    

    return X_new;
  }

}//namespace
