namespace ffd::heat_bath_mc{

  template<int max_order, int max_order_calc>
  std::array<Real, (1<<max_order)> conformal_map(std::array<Real, (1<<max_order)> const& P, std::vector<Real> const& map_matrix){

    BinaryInt const mat_size = ffd::core_math::sqrt_int(size(map_matrix));
    BinaryInt constexpr two_max_order = (1<<max_order);

    assert (mat_size*mat_size == (int)size(map_matrix));
    assert (max_order+1 == mat_size);

    auto map_mat = map_matrix;
      for (BinaryInt y=0; y<mat_size; ++y){
        for (BinaryInt x=0; x<=y; ++x){
        map_mat[x+y*mat_size] /= (Real)ffd::core_math::BinomialCoefficient(y, x);
      }
    }

    std::array<Real, two_max_order> R;

    for (BinaryInt V=1; V < two_max_order; ++V){
      BinaryInt const card_V_mat = __builtin_popcount(V) * mat_size;
      for (BinaryInt S = V; S > 0; S = ((S-1)&V)){
        BinaryInt const card_S = __builtin_popcount(S);
        R[V] += P[S] * map_mat[card_V_mat+card_S];
      }
    }
    R[0] = P[0];
    return R;
  }


}//namespace
