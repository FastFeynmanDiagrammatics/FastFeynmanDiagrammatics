namespace ffd::inverse_matrix_gauss::unit_test {

void one_by_one() {
  using std::abs;
  Real const eps_real = 100 * std::numeric_limits<Real>::epsilon();
  std::array<Real, 1> M;
  M.fill(3.1415);
  auto M_copy = M;
  auto M_copy2 = M;

  auto [inverse, det] = Inverse_and_Determinant(M);
  assert((abs(inverse[0] - 1 / M_copy[0]) < eps_real));
  auto [cofactor, det_c] = Cofactor_and_Determinant(M_copy2);
  assert((abs(cofactor[0] / det_c - 1 / M_copy[0]) < eps_real));

  auto det2 = ffd::determinant::Determinant(M_copy);
  assert((abs(det - det2) < eps_real));
  assert((abs(det_c - det2) < eps_real));

  // std::cerr<<1./3.1415<<std::endl;
  // std::cerr<<inverse[0]<<" "<<det<<"\n";
}

}  // namespace ffd::inverse_matrix_gauss::unit_test
