namespace ffd::sigmadet_hubbard {

//  ffd::user_space::Timer t_temp, t_1, t_2, t_3, t_4;

template <std::size_t n, class field, class poly_t>
void CDet_Xi_poly_r(
    array2d<field, n, n> const& __restrict__ G0,
    vector3d<poly_t>& __restrict__ Xi,
    array2d<std::uint8_t, n + 1, (1ul << n)> const& __restrict__ bm) {
  //    t_3.ini();
  static_assert(std::is_same_v<field, std::decay_t<decltype(Xi[0][0])>>);
  assert((size_x(Xi) == n));
  assert((size_y(Xi) == n));
  assert((size_z(Xi) == (1ul << n)));

  Xi.fill(field(0));
  array1d<poly_t, (1ul << n)> Z;
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
    //   auto const det2 = det * det;
    Z[S] = det;

    for (ulong j = 0; j < card_S; ++j) {
      auto const bm_j1_S = bm(j + 1, S);
      auto const j_S = j * card_S;
      for (ulong k = 0; k < card_S; ++k) {
        Xi(bm(k + 1, S), bm_j1_S, S) = C[j_S + k];
      }
    }
  }

  ulong const poly_n = Xi[0].coef.size() - 1;
  for (ulong S = 1; S < (1ul << n); ++S) {
    auto const bm_0_S = bm(0, S);
    auto const k_max = std::min(poly_n, n - bm_0_S);
    for (ulong V = ((S - 1) & S); V != 0; V = ((V - 1) & S)) {
      auto const bm_0_V = bm(0, V);
      auto const k = bm_0_S - bm_0_V;
      if (k <= k_max) {
        Z[S][k] += Z[V][0];
        for (ulong j = 0; j < bm_0_V; ++j) {
          auto const bm_j_V = bm(j + 1, V);
          for (ulong l = 0; l < bm_0_V; ++l) {
            auto const bm_l_V = bm(l + 1, V);
            Xi(bm_l_V, bm_j_V, S)[k] += Xi(bm_l_V, bm_j_V, V)[0];
          }  // for l in range(0, bm_0_S)
        }    // for j in range(0, bm_0_S)
      }      // if k
    }        // for V
    if (bm_0_S <= k_max) {
      Z[S][bm_0_S] += Z[0][0];
    }
  }  // for S in range(1, (1ul<<n))

  for (ulong S = 1; S < (1ul << n); ++S) {
    auto const bm_0_S = bm(0, S);
    bool const is_even_S = ((bm_0_S & 1) == 0);
    for (ulong j = 0; j < bm_0_S; ++j) {
      auto const bm_j1_S = bm(j + 1, S);
      for (ulong k = 0; k < bm_0_S; ++k) {
        Xi(bm(k + 1, S), bm_j1_S, S) *= Z[S];
        if (is_even_S) {
          Xi(bm(k + 1, S), bm_j1_S, S) = -Xi(bm(k + 1, S), bm_j1_S, S);
        }
      }
    }
    Z[S] *= Z[S];
    if (!is_even_S)
      Z[S] = -Z[S];
  }  // for S in range(0, (1ul<<n))// compute connected part

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
}  // routine CDet_Xi_poly

}  // namespace ffd::sigmadet_hubbard
