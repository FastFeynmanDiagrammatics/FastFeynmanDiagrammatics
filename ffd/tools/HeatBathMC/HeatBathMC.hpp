namespace ffd::heat_bath_mc{
  
  template<int max_order,
	   int target_order,
	   bool is_translation_invariant,
	   int max_calculated_order,
	   typename poly_t,
	   typename seed_t,
	   typename vector_coord_t,
	   typename coord_t>

  std::pair<std::array<Real, (1<<max_order)>,
    std::array<ffd::type_traits::element_of_container_t<vector_coord_t>, max_order>>

  HeatBathMC(poly_t const& P,
	     seed_t const& Phi,
	     vector_coord_t const& X,
	     coord_t const& O){
    static_assert(
		  (is_translation_invariant &&
		   std::is_same_v<typename std::decay_t<coord_t>, default_tmp_t>)
		  || (!is_translation_invariant &&
		      !std::is_same_v<typename std::decay_t<coord_t>, default_tmp_t>)
		  );
    assert(( size(X) == max_order ));
    using namespace ffd::type_traits;
    if constexpr(is_std_array_v<std::decay_t<vector_coord_t>>){
         static_assert(constexpr_size_v<std::decay_t<vector_coord_t>> == max_order);
}
    
    
    int constexpr eff_max_order = max_order + !is_translation_invariant;
    using element_t = typename ffd::type_traits::element_of_container_t<vector_coord_t>;
    std::array<element_t, eff_max_order> X_extended;
    for(int j=0; j<max_order; ++j){
      X_extended[j] = X[j];
    }
    if constexpr(!is_translation_invariant){
	static_assert(std::is_same_v<coord_t, element_t>);
	X_extended[max_order] = O;
      }
    

    auto const seeds = compute_seeds<eff_max_order>(Phi, X_extended);
    

    auto const P_cond = compute_conditional_probabilities<max_order,
							  is_translation_invariant>(seeds);
    
    
    auto const fluxes = compute_fluxes<max_order, max_calculated_order>(P, P_cond);
   

    auto const [times, chosen_set] =
      heat_bath_fluxes<max_order, target_order, max_calculated_order>(fluxes);
        

    assert(( __builtin_popcount(chosen_set) == target_order ));

    
    auto const X_new = create_new_vertices<max_order,
					   is_translation_invariant>(Phi,
								     X_extended,
								     chosen_set);
        
    return std::make_pair(times, X_new);
  }

}//namespace
