namespace ffd::nilpotent_polynomial_s::unit_test {
template <auto n, mode m = safe>
void test_multiply_divide() {
  P<Real, n, m> P0, P1;
  P0.precision *= 1;
  for (ulong j = 0; j < (1ul << n); ++j) {
    P0[j] = 2 * std::cos(10 * j);
    P1[j] = (j + 1) * std::sin(13 * j + 1);
  }  // for j in range(0, (1ul<<n))
  auto P2 = (P0 * P1) / P1 - P0 * 2. + P0;
  for (ulong j = 0; j < (1ul << n); ++j) {
    std::cerr << j << " " << P2[j] << "\n";
  }  // for j in range(0, (1ul<<n))
  std::cerr << "\n\n";
  P2 = (P0 * P1) / P1 / P0;
  for (ulong j = 0; j < (1ul << n); ++j) {
    std::cerr << j << " " << P2[j] << "\n";
  }  // for j in range(0, (1ul<<n))
  std::cerr << "\n\n";
  std::cerr << "\n\n";
}
}  // namespace ffd::nilpotent_polynomial_s::unit_test
