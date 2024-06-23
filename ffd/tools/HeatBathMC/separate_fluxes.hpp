namespace ffd::heat_bath_mc{

  template<int target_order,
	   int max_order,
	   typename flux_t>

  auto
  separate_fluxes(flux_t const& fluxes){
    BinaryInt constexpr two_max_order = (1<<max_order);
    flux_t fluxes_u = fluxes;
    flux_t fluxes_s = fluxes;

    
    for(BinaryInt S=0; S < two_max_order; ++S){
      auto const card_S = __builtin_popcount(S);
      if( card_S != target_order ){
	fluxes_u[S] = 0.;
      }else{
	fluxes_s[S] = 0.;
      }
    }
    
    
    return std::make_pair(fluxes_u, fluxes_s);
  }

}//namespace
