namespace ffd::heat_bath_mc{

  template<int max_order, int max_order_calc>
  std::pair<std::array<Real, (1<<max_order)>, BinaryInt>
  heat_bath_fluxes_map(std::array<Real, (1<<max_order)> const& fluxes, std::vector<Real> const& map_matrix){

    static const std::array<Real, max_order_calc+1> one_bin_coef = compute_bin_coef<max_order, max_order_calc>();
    int constexpr minus_order = (max_order_calc == max_order);

    BinaryInt const mat_size = ffd::core_math::sqrt_int(size(map_matrix));
    BinaryInt constexpr two_max_order = (1<<max_order);
    assert (mat_size*mat_size == (int)size(map_matrix));
    assert (max_order+1 == mat_size);

    std::array<Real, max_order+1> fac;
    fac[0]=1.;
    for (BinaryInt j=0; j<=max_order; ++j){
      fac[j+1] = (j+1) * fac[j];
    }

    auto map_mat = map_matrix;
      for (BinaryInt y=0; y<mat_size; ++y){
        for (BinaryInt x=0; x<=y; ++x){
          map_mat[x+y*mat_size] *= fac[y-x];
      }
    }

    std::array<Real, two_max_order> mapped_fluxes;
    mapped_fluxes.fill(0.);
    mapped_fluxes[0] = fluxes[0];
    for (BinaryInt V = 1; V < two_max_order; ++V){
      BinaryInt const card_V_mat = __builtin_popcount(V) * mat_size;
      for (BinaryInt S = V; S > 0; S = ((S-1)&V)){
        BinaryInt const card_S = __builtin_popcount(S);
        mapped_fluxes[V] += fluxes[S] * map_mat[card_V_mat+card_S];
      }
    }

    std::array<Real, (1<<max_order)> normed_fluxes;
    normed_fluxes.fill(0.);
    Real w_target = 0.;
    for(BinaryInt S=0; S<(1<<max_order_calc)-minus_order; ++S){
        normed_fluxes[S] = std::abs(mapped_fluxes[S])*one_bin_coef[__builtin_popcount(S)];
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
	      times[S] = mapped_fluxes[S]*one_over_w_target;
      }
    }

    return std::make_pair(times, chosen_S);
  }


  template<int max_order, int max_order_calc>
  std::pair<std::array<Real, (1<<max_order)>, BinaryInt>
  heat_bath_fluxes_map_post(std::array<Real, (1<<max_order)> const& fluxes, std::vector<Real> const& map_matrix){

    static const std::array<Real, max_order_calc+1> one_bin_coef = compute_bin_coef<max_order, max_order_calc>();
    int constexpr minus_order = (max_order_calc == max_order);

    BinaryInt const mat_size = ffd::core_math::sqrt_int(size(map_matrix));
    BinaryInt constexpr two_max_order = (1<<max_order);
    assert (mat_size*mat_size == (int)size(map_matrix));
    assert (max_order+1 == mat_size);

    std::array<Real, max_order+1> fac;
    fac[0]=1.;
    for (BinaryInt j=0; j<=max_order; ++j){
      fac[j+1] = (j+1) * fac[j];
    }

    auto map_mat = map_matrix;
      for (BinaryInt y=0; y<mat_size; ++y){
        for (BinaryInt x=0; x<=y; ++x){
          map_mat[x+y*mat_size] *= fac[y-x];
      }
    }

    std::array<Real, (1<<max_order)> normed_fluxes;
    normed_fluxes.fill(0.);
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

    std::array<Real, two_max_order> mapped_times;
    mapped_times.fill(0.);
    mapped_times[0] = times[0];
    for (BinaryInt V = 1; V < two_max_order; ++V){
      BinaryInt const card_V_mat = __builtin_popcount(V) * mat_size;
      for (BinaryInt S = V; S > 0; S = ((S-1)&V)){
        BinaryInt const card_S = __builtin_popcount(S);
        mapped_times[V] += times[S] * map_mat[card_V_mat+card_S];
      }
    }

    return std::make_pair(mapped_times, chosen_S);
  }

}//namespace
