namespace ffd::heat_bath_mc{

  template<int max_order, int max_order_calc>
  std::pair<std::array<Real, (1<<max_order)>, BinaryInt>
  heat_bath_fluxes_sum(std::array<Real, (1<<max_order)> const& fluxes, Real const U=1.){

    static const std::array<Real, max_order_calc+1> one_bin_coef = compute_bin_coef<max_order, max_order_calc>();
    int constexpr minus_order = (max_order_calc == max_order);

    std::array<Real, max_order+1> powU;
    powU[0]=1.;
    for (BinaryInt j=0; j<max_order; ++j){
      powU[j+1] = U * powU[j];
    }

    std::array<Real, (1<<max_order)> summed_fluxes;
    summed_fluxes.fill(0.);
    for (BinaryInt S=0; S<(1<<max_order); ++S){
      for (BinaryInt j=1; j<S; j*=2){
        summed_fluxes[S] += ((j&S)==j) * summed_fluxes[S-j];
      }
      summed_fluxes[S] += fluxes[S] * powU[__builtin_popcount(S)];
    }
    
    std::array<Real, (1<<max_order)> normed_fluxes;
    normed_fluxes.fill(0.);
    Real w_target = 0.;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
        normed_fluxes[S] = std::abs(summed_fluxes[S])*one_bin_coef[__builtin_popcount(S)];
      w_target += normed_fluxes[S];
    }

    Real const pw_target = w_target*ffd::user_space::Proba();
    Real pw_cumul = 0.;
    BinaryInt chosen_S;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
      pw_cumul += normed_fluxes[S];
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
	      times[S] = summed_fluxes[S]*one_over_w_target;
      }
    }

    return std::make_pair(times, chosen_S);
  }


  template<int max_order, int max_order_calc>
  std::pair<std::array<Real, (1<<max_order)>, BinaryInt>
  heat_bath_fluxes_sum_post(std::array<Real, (1<<max_order)> const& fluxes, Real const U=1.){

    static const std::array<Real, max_order_calc+1> one_bin_coef = compute_bin_coef<max_order, max_order_calc>();
    int constexpr minus_order = (max_order_calc == max_order);

    std::array<Real, max_order+1> powU;
    powU[0]=1.;
    for (BinaryInt j=0; j<max_order; ++j){
      powU[j+1] = U * powU[j];
    }

    std::array<Real, (1<<max_order)> normed_fluxes;

    Real w_target = 0.;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
      normed_fluxes[S] = std::abs(fluxes[S])*one_bin_coef[__builtin_popcount(S)];
      w_target += normed_fluxes[S];
    }

    Real const pw_target = w_target*ffd::user_space::Proba();
    Real pw_cumul = 0.;
    BinaryInt chosen_S;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
      pw_cumul += normed_fluxes[S];
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

    std::array<Real, (1<<max_order)> summed_times;
    summed_times.fill(0.);
    for (BinaryInt S=0; S<(1<<max_order); ++S){
      for (BinaryInt j=1; j<S; j*=2){
        summed_times[S] += ((j&S)==j) * summed_times[S-j];
      }
      summed_times[S] += times[S]*powU[__builtin_popcount(S)];
    }

    return std::make_pair(summed_times, chosen_S);
  }

}//namespace
