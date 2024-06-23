namespace ffd::inverse_matrix_gauss::unit_test {

template <auto n, auto poly_n, bool non_zero_diagonal>
auto taylor_poly() {
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
      if (j == k)
        M_a(j, k)[1] = Real(1);
    }  // for k in range(0, n)
  }  // for j in range(0, n)

  auto M_cp = M;
  auto M_a_cp = M_a;

  auto [I, det] = Inverse_and_Determinant(M);
  auto [I_a, det_a] = Inverse_and_Determinant(M_a);

  for (ulong j = 0; j < n; ++j) {
    for (ulong k = 0; k < n; ++k) {
      auto prod = M_cp(j, 0) * I(0, k);
      auto prod_a = M_a_cp(j, 0) * I_a(0, k);
      for (ulong l = 1; l < n; ++l) {
        prod += M_cp(j, l) * I(l, k);
        prod_a += M_a_cp(j, l) * I_a(l, k);
      }
      //      std::cerr << (j == k) << " == ";
      for (ulong p = 0; p < poly_n + 1; ++p) {
        //        std::cerr << prod_a[p] << " ";
      }  // for p in range(0, poly_n+1)prod_a
      //      std::cerr << "\n";
    }
  }
  //  std::cerr << "det = " << det << " " << det_a[0] << "\n";
  for (ulong j = 0; j < n; ++j) {
    for (ulong k = 0; k < n; ++k) {
      //  std::cerr << j << " " << k << " =" << I(j, k) * det << " "
      //      << (I_a(j, k) * det_a)[0] << "\n";
      for (ulong p = 1; p < poly_n + 1; ++p) {
        //        std::cerr << j << " " << k << " " << p << " " << (I_a(j, k) *
        //        det_a)[p]
        //          << "\n";

      }  // for p in range(1, poly_n+1)
    }  // for k in range(0, n)
  }  // for j in range(0, n)
}

}  // namespace ffd::inverse_matrix_gauss::unit_test
