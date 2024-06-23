namespace ffd::heat_bath_mc{

  template<int max_order,
	   int max_order_calc>

  std::pair<std::array<Real, (1<<max_order)>, BinaryInt>

  heat_bath_fluxes_all(std::array<Real, (1<<max_order)> const& fluxes){

    static const std::array<Real, max_order_calc+1> one_bin_coef = compute_bin_coef<max_order, max_order_calc>();

    int constexpr minus_order = (max_order_calc == max_order);

    Real w_target = 0.;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
      w_target += std::abs(fluxes[S])*one_bin_coef[__builtin_popcount(S)];
    }


    Real const pw_target = w_target*ffd::user_space::Proba();
    Real pw_cumul = 0.;
    BinaryInt chosen_S;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
      pw_cumul += std::abs(fluxes[S])*one_bin_coef[__builtin_popcount(S)];
      if(pw_cumul > pw_target){
	      chosen_S = S;
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
