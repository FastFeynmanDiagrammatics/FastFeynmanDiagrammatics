namespace ffd::pade::unit_test {

void simple_rational() {
  std::vector<Real> coef(20);
  Real pow_R = 1, A = 3.14, R = -2;
  for (ulong j = 0; j < size(coef); ++j, pow_R *= R) {
    coef[j] = A * pow_R;
  }  // for j in range(0, size(coef))]
  int p_max = 2;
  int q_max = 3;
  auto [p, q] = eval(p_max, q_max, coef);
  for (ulong j = 0; j < size(p); ++j) {
    std::cerr << "p " << j << " " << p[j] << "\n";
  }  // for j in range(0, size(p))
  for (ulong j = 0; j < size(q); ++j) {
    std::cerr << "q " << j + 1 << " " << q[j] << "\n";
  }  // for j in range(0, size(p))
  for (ulong j = 0; j < p_max + q_max + 1; ++j) {
    std::cerr << "c " << j << " " << coef[j] << "\n";
  }  // for j in range(0, size(p))
}

}  // namespace ffd::pade::unit_test
