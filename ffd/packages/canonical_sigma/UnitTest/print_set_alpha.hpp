namespace ffd::canonical_sigma::unit_test {

void print_set_alpha() {
  int constexpr n = 5;
  ffd::array1d<Real, (1ul << n)> PM;
  ffd::vector2d<Real> pow_alpha;
  for (ulong S = 0; S < (1ul << n); ++S) {
    PM[S] = std::cos(S * 70);
  }  // for S in range(0, (1ul<<n))
  PM[0] = 1;
  for (ulong S = 0; S < (1ul << n); ++S) {
    if (__builtin_popcount(S) == 1)
      PM[S] = 0;
  }  // for S in range(0, (1ul<<n))
  auto A0 = set_alpha_r<n, n>(PM, pow_alpha);
  auto A1 = set_alpha_r<n, n - 1>(PM, pow_alpha);
  auto A2 = set_alpha_r<n, n - 2>(PM, pow_alpha);
  auto A3 = set_alpha_r<n, n - 3>(PM, pow_alpha);
  auto A4 = set_alpha_r<n, n - 3>(PM, pow_alpha);
  auto A5 = set_alpha_r<n, n - 3>(PM, pow_alpha);
  for (ulong S = 0; S < (1ul << n); ++S) {
    for (ulong j = 0; j < n + 1; ++j) {
      std::cerr << S << " " << j << " " << A0[S][j] << ", err= ";
      if (j < A1[S].coef.size())
        std::cerr << A0[S][j] - A1[S][j] << " ";
      if (j < A1[S].coef.size())
        std::cerr << A0[S][j] - A1[S][j] << " ";
      if (j < A2[S].coef.size())
        std::cerr << A0[S][j] - A2[S][j] << " ";
      if (j < A3[S].coef.size())
        std::cerr << A0[S][j] - A3[S][j] << " ";
      if (j < A4[S].coef.size())
        std::cerr << A0[S][j] - A4[S][j] << " ";
      if (j < A5[S].coef.size())
        std::cerr << A0[S][j] - A5[S][j] << " ";
      std::cerr << ", order = " << j + __builtin_popcount(S) << "\n";
    }  // for j in range(0, half_n+1)
  }    // for S in range(0, (1ul<<n))
}

}  // namespace ffd::canonical_sigma::unit_test
