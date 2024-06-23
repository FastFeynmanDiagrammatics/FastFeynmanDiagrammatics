namespace ffd::inverse_matrix_gauss::unit_test {

template <int n>
void n_by_n() {
  using std::abs;
  long double const eps_real =
      100000 * std::numeric_limits<long double>::epsilon();
  std::array<long double, n * n> M;
  for (int j = 0; j < n * n; ++j) {
    M[j] = 2 * ffd::user_space::Proba() - 1.;
  }
  auto M_copy = M;
  auto M_copy2 = M;

  auto [inverse, det] = Inverse_and_Determinant(M);
  auto [cofactor, det_c] = Cofactor_and_Determinant(M_copy2);
  decltype(M) identity, identity_c;
  identity.fill(0.);
  identity_c.fill(0.);
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      for (int l = 0; l < n; ++l) {
        identity[j * n + k] += inverse[j * n + l] * M_copy[l * n + k];
        identity_c[j * n + k] += cofactor[j * n + l] * M_copy[l * n + k];
      }
      // std::cerr<<std::setprecision(1)<<identity[j*n+k]<<" ";
      // if( k == n-1 ){
      //   std::cerr<<"\n";
      // }

      assert((abs(identity[j * n + k] - (j == k)) < eps_real));
      assert((abs(identity_c[j * n + k] / det_c - (j == k)) < eps_real));
    }
  }

  auto det2 = ffd::determinant::Determinant(M_copy);
  // std::cerr<<std::setprecision(20)<<"det = "<<det2<<" "<<det<<std::endl;
  assert((abs(det - det2) < eps_real));
  assert((abs(det_c - det2) < eps_real));
}

}  // namespace ffd::inverse_matrix_gauss::unit_test
