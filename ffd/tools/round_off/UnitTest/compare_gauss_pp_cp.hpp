namespace ffd::round_off::unit_test {

template <int n>
void compare_gauss_pp_cp() {
  using std::abs;
  using real_t = R<Real>;
  ffd::array2d<real_t, n, n> M;
  M.fill(0);
  bool use_exp = false;
  for (ulong j = 0; j < 1ul << 20; ++j) {
    auto rand_num = [] { return 2 * ffd::user_space::Proba() - 1.; };
    for (ulong j = 0; j < n; ++j) {
      for (ulong k = 0; k < n; ++k) {
        M(j, k) = rand_num() *
                  (use_exp ? std::pow(10, -15 * (1 + rand_num()) * .5) : 1);
      }  // for k in range(0, n)
    }    // for j in range(0, n)
    auto M_cp = M;

    auto det = ffd::determinant::Determinant(M);
    auto det_cp = ffd::determinant::Gauss_CP(M_cp);
    auto [det_val, det_err] = det;
    //    std::cerr << std::setprecision(15) << (det - det_cp).val << " " <<
    //    det_err
    //        << " " << det_cp.val << "\n";
    //    std::cerr << det_val << " +/- " << det_err << " here\n";
    assert((abs(det - det_cp) <= det_err));
  }  // for j in range(0, 1ul<<10)
}
}  // namespace ffd::round_off::unit_test
