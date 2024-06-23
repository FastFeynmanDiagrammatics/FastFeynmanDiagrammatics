namespace ffd::heat_bath_mc{

  template<int max_order,
	   int target_order,
	   int max_order_calc>

  std::pair<std::array<Real, (1<<max_order)>, BinaryInt>

  heat_bath_fluxes(std::array<Real, (1<<max_order)> const& fluxes){
    int constexpr num_conf_target =
      ffd::core_math::BinomialCoefficient<max_order, target_order>();
    std::array<Real, num_conf_target> p_target;
    std::array<BinaryInt, num_conf_target> S_target;
    
    
    int target_counter = 0;
    for(BinaryInt S=0; S<(1<<max_order); ++S){
      if(__builtin_popcount(S) == target_order){
	p_target[target_counter] = std::abs( fluxes[S] );
	S_target[target_counter] = S;
	++target_counter;
      }
    }
    
    
    Real w_target = 0.;
    for(int j=0; j<num_conf_target; ++j){
      w_target += p_target[j];
    }

    
    Real const pw_target = w_target*ffd::user_space::Proba();
    Real pw_cumul = 0.;
    BinaryInt chosen_S;
    for(int j=0; j<num_conf_target; ++j){
      pw_cumul += p_target[j];
      if(pw_cumul > pw_target){
	chosen_S = S_target[j];
	break;
      }
    }
    

    Real const one_over_w_target = 1./w_target;
    std::array<Real, (1<<max_order)> times;
    times.fill(0.);
    for(BinaryInt S=0; S<(1<<max_order); ++S){
      if(__builtin_popcount(S) <= max_order_calc){
	times[S] = fluxes[S]*one_over_w_target;
      }
    }
    
    
    return std::make_pair(times, chosen_S);
  }



}//namespace
