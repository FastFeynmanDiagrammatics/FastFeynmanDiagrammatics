namespace ffd::heat_bath_mc{

  template<int eff_order,
	   typename seed_t,
	   typename vector_coord_t>

  std::array<Real, eff_order*eff_order>

  compute_seeds(seed_t const& Phi,
		vector_coord_t const& X){
    using namespace ffd::type_traits;
    if constexpr(is_std_array_v<std::decay_t<vector_coord_t>>){
        static_assert(constexpr_size_v<std::decay_t<vector_coord_t>> == eff_order);
      }
    assert(( unsigned(eff_order) == size(X) ));
    std::array<Real, eff_order*eff_order> seeds;
    
    
    for(int j=0; j<eff_order; ++j){
      for(int k=0; k<eff_order; ++k){
	if( j != k ){
	  seeds[j*eff_order+k] = Phi(X[j], X[k]);
	}else{
	  seeds[j*eff_order+k] = 0.;
	}
      }
    }
    
    
    return seeds;
  }
  
  
}//namespace
