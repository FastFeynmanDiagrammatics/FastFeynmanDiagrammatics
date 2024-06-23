namespace ffd::sigmadet_hubbard {

//  ffd::user_space::Timer t_temp, t_1, t_2, t_3, t_4;

template <std::size_t n, class field = Real>
void CDet_Xi_r(
    array2d<field, n, n> const& __restrict__ G0,
    vector3d<field>& __restrict__ Xi,
    array2d<std::uint8_t, n + 1, (1ul << n)> const& __restrict__ bm) {
  //    t_3.ini();
  assert((size_x(Xi) == n));
  assert((size_y(Xi) == n));
  assert((size_z(Xi) == (1ul << n)));

  Xi.fill(field(0));
  array1d<field, (1ul << n)> Z;
  Z.fill(field(0));
  Z[0] = field(1);
  for (BinaryInt S = 1; S < (1ul << n); ++S) {
    ulong const card_S = bm(0, S);
    bool const is_even_S = ((card_S & 1) == 0);

    vector2d<field> g0(card_S, card_S);
    //  t_2.ini();
    for (ulong j = 0; j < card_S; ++j) {
      auto const bm_j1_S = bm(j + 1, S);
      for (ulong k = 0; k < card_S; ++k) {
        g0(k, j) = G0(bm(k + 1, S), bm_j1_S);
      }
    }
    // t_2.fin();
    //      t_temp.ini();
    auto const [C, det] =
        ffd::inverse_matrix_gauss::Cofactor_and_Determinant(std::move(g0));
    //      t_temp.fin();
    auto const det2 = det * det;
    Z[S] = is_even_S ? det2 : -det2;

    for (ulong j = 0; j < card_S; ++j) {
      auto const bm_j1_S = bm(j + 1, S);
      for (ulong k = 0; k < card_S; ++k) {
        auto const xi_temp = C[j * card_S + k] * det;
        Xi(bm(k + 1, S), bm_j1_S, S) = is_even_S ? -xi_temp : xi_temp;
      }
    }
  }
  // compute connected part
  //    t_1.ini();
  ulong constexpr n2 = n * n;
  for (ulong V = 1; V < (1ul << n); ++V) {
    auto const V_n = n2 * V;
    for (ulong S = ((V - 1) & V); S != 0; S = ((S - 1) & V)) {
      auto const S_n = n2 * S;
      auto const Z_V_S = Z[V - S];
      ulong const bm_0_S = bm(0, S);
      for (ulong j = 0; j < bm_0_S; ++j) {
        ulong const bmj_n = bm(j + 1, S) * n;
        ulong const V_j_n = V_n + bmj_n;
        ulong const S_j_n = S_n + bmj_n;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (ulong k = 0; k < bm_0_S; ++k) {
          auto const bm_k1_S = bm(k + 1, S);
          Xi[V_j_n + bm_k1_S] -= Xi[S_j_n + bm_k1_S] * Z_V_S;
        }
      }  // for j in range(0, bm(0, S))
    }
  }
  // t_1.fin();
  // t_3.fin();
  return;
}

}  // namespace ffd::sigmadet_hubbard
