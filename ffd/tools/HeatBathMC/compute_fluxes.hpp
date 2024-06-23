namespace ffd::heat_bath_mc{

  template<int max_order,
	   int max_order_calc,
	   typename poly_t,
	   class array_t = int>
  
  std::array<Real, (1<<max_order)>
  
  compute_fluxes(poly_t const& P,
		 std::array<Real, (1<<max_order)> const& P_cond,
		 array_t const& lambda_S = 1){
    std::array<Real, max_order_calc+1> lambdas;
    lambdas.fill(1.);
    if constexpr(!std::is_same_v<array_t, int>){
	for(uint j=0; j<max_order+1; ++j){
	  lambdas[j] = lambda_S[j];
	}
      }
    using std::abs;
    BinaryInt constexpr two_max_order = (1<<max_order);
    std::array<Real, two_max_order> fluxes;
    fluxes.fill(0.);


    for(BinaryInt S=0; S < two_max_order; ++S){
      auto pop = __builtin_popcount(S);
      if(pop <= max_order_calc){
	fluxes[S] = P[S]*P_cond[two_max_order-1-S]*lambdas[pop];
      }
    }
    

    return fluxes;
  }

}//namespace
