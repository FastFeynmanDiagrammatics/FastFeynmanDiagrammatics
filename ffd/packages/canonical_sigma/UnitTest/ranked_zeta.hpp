namespace ffd::canonical_sigma::unit_test {

void ranked_zeta() {
  int constexpr n = 5;
  ffd::array1d<Real, (1ul << n)> f;
  f.fill(0);
  for (ulong S = 0; S < (1ul << n); ++S) {
    f[S] = std::cos(40 * S);
  }  // for S in range(0, (1ul<<n))
  ffd::array2d<Real, (1ul << n), (n + 1)> f2;
  for (ulong S = 0; S < (1ul << n); ++S) {
    f2[S] = f[S];
  }  // for S in range(0, (1ul<<n))

  auto rkt = ranked_zeta_transform<n>(f);
  // ranked_zeta_transform_index_r<n>(f2, n);

  for (ulong S = 0; S < (1ul << n); ++S) {
    std::cerr << " " << S << " " << f(S) << "\n";
  }

  for (ulong S = 0; S < (1ul << n); ++S) {
    for (ulong j = 0; j < n + 1; ++j) {
      std::cerr << S << " " << j << " " << rkt(S, j) << " " << f2(S, j) << "\n";
    }  // for j in range(0, n+1)
  }  // for S in range(0, (1ul<<n))
}

}  // namespace ffd::canonical_sigma::unit_test
