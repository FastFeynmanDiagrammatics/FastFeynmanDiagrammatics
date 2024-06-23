namespace ffd::heat_bath_mc{

  
  template<int max_order,
	   bool is_translation_invariant,
	   typename seeds_t>

  auto
  single_pseudo_potential(seeds_t const& seeds,
			  BinaryInt S,
			  int vertex){
    assert((S&(1<<vertex))==0);
    std::vector<BinaryInt> const vector_S = ffd::set_theory::
      VectorOfBinaryDigitsOf(S);
    BinaryInt const card_S = size(vector_S);
    BinaryInt const max_order_eff = max_order + !is_translation_invariant;
    
    
    Real pseudo_pot = 0.;
    for(int j=0; j<card_S; ++j){
      pseudo_pot += seeds[vertex*max_order_eff +vector_S[j]];
    }
    if constexpr(!is_translation_invariant){
	pseudo_pot += seeds[vertex*max_order_eff +max_order];
      }else if(S==0){
      return 1.;
    }
    
    
    return pseudo_pot/(card_S+!is_translation_invariant);
  }

}//namespace
