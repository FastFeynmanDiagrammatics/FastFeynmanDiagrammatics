namespace ffd::inverse_matrix_gauss::unit_test {

template <auto n, auto poly_n, bool non_zero_diagonal>
auto taylor_poly_rzt() {
  using poly_t = ffd::taylor_polynomial::P<Real, poly_n>;

  ffd::array2d<Real, n, n> M;
  for (ulong j = 0; j < n; ++j) {
    for (ulong k = 0; k < n; ++k) {
      M(j, k) =
          (2 * ffd::user_space::Proba() - 1) * (j != k || non_zero_diagonal);
    }
  }  // for j in sizerange(M)

  ffd::array2d<poly_t, n, n> M_a;
  for (ulong j = 0; j < n; ++j) {
    for (ulong k = 0; k < n; ++k) {
      M_a(j, k) = M(j, k);
    }  // for k in range(0, n)
    M_a(j, j)[1] = Real(1);
  }  // for j in range(0, n)

  ffd::array1d<Real, (1ul << n)> Z;
  Z.fill(0);
  Z[0] = 1;
  ffd::array1d<Real, n, n, (1ul << n)> A;
  A.fill(0);
  auto const bm = for (ulong S = 1; S < (1ul << n); ++S) {
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
    auto const [I, det] =
        ffd::inverse_matrix_gauss::Inverse_and_Determinant(std::move(g0));
    //      t_temp.fin();
    auto const det2 = det * det;
    Z[S] = is_even_S ? det2 : -det2;

    for (ulong j = 0; j < card_S; ++j) {
      auto const bm_j1_S = bm(j + 1, S);
      for (ulong k = 0; k < card_S; ++k) {
        auto const xi_temp = I[j * card_S + k] * det2;
        Xi(bm(k + 1, S), bm_j1_S, S) = is_even_S ? -xi_temp : xi_temp;
      }
    }
  }
}

}  // namespace ffd::inverse_matrix_gauss::unit_test
