namespace ffd::determinant::unit_test {

template <int n>
void compare_gauss_quad() {
  using real_t = double;
  ffd::array2d<real_t, n, n> M;
  M.fill(0);
  bool use_exp = false;

  auto rand_num = [] { return 2 * ffd::user_space::Proba() - 1.; };
  for (ulong j = 0; j < n; ++j) {
    for (ulong k = 0; k < n; ++k) {
      M(j, k) = rand_num() *
                (use_exp ? std::pow(10, -15 * (1 + rand_num()) * .5) : 1);
    }  // for k in range(0, n)
  }  // for j in range(0, n)
  auto M_cp = M;
  auto M_cp2 = M;

  auto det = Determinant(M);
  auto det_cp = Gauss_CP(M_cp);
  auto det_cp2 = Gauss_CP_quad(M_cp2);
  double det_cp2_double = det_cp2;
  //  std::cerr << (det - det_cp) / det_cp << " " << det_cp << " here\n";
  //  std::cerr << std::setprecision(16) << det << "\n"
  //        << det_cp << "\n"
  //          << det_cp2_double << "\n";
  //  printf("%Q\nhere\n", det_cp2);
  assert((std::abs((det - det_cp) / det_cp) < 1e-13));
}
}  // namespace ffd::determinant::unit_test
