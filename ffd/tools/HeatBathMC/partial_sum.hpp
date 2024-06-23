namespace ffd::heat_bath_mc{

  template<int max_order, int max_order_calc>
  std::array<Real, (1<<max_order)> partial_sum(std::array<Real, (1<<max_order)> const& P, Real const U=1.){

   //  static const std::array<Real, max_order_calc+1> one_bin_coef =
   //    compute_bin_coef<max_order, max_order_calc>();
   // for (uint i=0; i<=max_order_calc; ++i){
   // }

    BinaryInt constexpr two_max_order = (1<<max_order);

    std::array<Real, two_max_order> R;
    R.fill(0.);

    std::array<Real, max_order+1> powU;
    powU[0]=1.;
    for (BinaryInt j=0; j<max_order; ++j){
      powU[j+1] = U * powU[j];
    }

    for (BinaryInt S=1; S<two_max_order; ++S){
      BinaryInt const card_S = __builtin_popcount(S);
      auto temp = P[S];
      // auto temp = P[S] * one_bin_coef[card_S];
      for (BinaryInt j=1; j<S; j*=2){
        R[S] += ((j&S)==j) * R[S-j];
      }
      R[S] /= double(card_S);
      R[S] += powU[card_S]*temp;
    }
    R[0] = P[0];

    // exit(1);
    return R;
  }

}//namespace
