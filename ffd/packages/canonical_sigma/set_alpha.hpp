namespace ffd::canonical_sigma {

template <std::size_t n, std::size_t half_n, class vector_t, class field_d>

auto set_alpha_r(vector_t const& PM,
                 ffd::vector2d<field_d>& __restrict__ pow_alpha) {
  //  std::size_t constexpr half_n = n;  // (n - (n & 1)) / 2;
  using tpoly_t = ffd::truncated_polynomial::P<half_n, Real>;
  //  auto const PM = ffd::principal_minors::PrincipalMinorsBFS(G0);
  auto const rzt = ranked_zeta_transform<n>(PM);
  ffd::vector1d<tpoly_t> PM_a(1ul << n);
  for (ulong S = 0; S < (1ul << n); ++S) {
    for (ulong j = 0; j <= std::min<int>(half_n, __builtin_popcount(S)); ++j) {
      PM_a[S][j] = rzt(S, __builtin_popcount(S) - j);
    }  // for j in range(0, half_n)
  }    // for S in range(0, (1ul<<n))
  ffd::nilpotent_polynomial::NilpotentPolynomial<tpoly_t> A(n), Z(n);
  for (ulong S = 0; S < (1ul << n); ++S) {
    Z[S] = 0.;
    A[S] = 0.;
    if ((__builtin_popcount(S) & 1) == 0) {
      Z[S] = PM_a[S] * PM_a[S];
      for (ulong j = 0; j < n; ++j) {
        if (((S >> j) & 1) != 0) {
          A[S] -= PM_a[S - (1ul << j)];
        }
      }  // for j in range(0, order)
    } else {
      Z[S] = -PM_a[S] * PM_a[S];
      for (ulong j = 0; j < n; ++j) {
        if (((S >> j) & 1) != 0) {
          A[S] += PM_a[S - (1ul << j)];
        }
      }  // for j in range(0, order)
    }    // if else
    A[S] = A[S] * PM_a[S];
  }  // for S in range(0, (1ul<<n))
  A /= Z;
  for (uint S = 1; S < (1ul << n); S = ffd::set_theory::next_permutation(S)) {
    A[S] = 0.;
  }  // for S
  A /= Real(n);

  for (ulong k = 2; k < n; ++k) {
    for (uint S = (1ul << k) - 1; S < (1ul << n);
         S = ffd::set_theory::next_permutation(S)) {
      ;
    }  // for S
  }    // for k

  return A;
}
}  // namespace ffd::canonical_sigma
