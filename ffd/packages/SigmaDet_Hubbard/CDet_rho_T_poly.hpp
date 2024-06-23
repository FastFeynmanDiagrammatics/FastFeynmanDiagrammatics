namespace ffd::sigmadet_hubbard {

template <std::size_t n, class field, class poly_t>
void CDet_rho_T_poly_r(
    array2d<field, n, n> const& __restrict__ G0,
    vector3d<poly_t> const& __restrict__ Xi,
    vector3d<poly_t>& __restrict__ rho_T,
    array2d<std::uint8_t, n + 1, (1ul << n)> const& __restrict__ bm) {
  assert((size_x(Xi) == n));
  assert((size_y(Xi) == n));
  assert((size_z(Xi) == (1ul << n)));
  assert((size_x(rho_T) == n));
  assert((size_y(rho_T) == n));
  assert((size_z(rho_T) == (1ul << n)));

  rho_T.fill(field(0));
  ulong constexpr n2 = n * n;
  for (ulong S = 1; S < (1ul << n); ++S) {
    ulong const bm_0_S = bm(0, S);
    ulong const S_n = n2 * S;
    for (ulong k = 1; k <= bm_0_S; ++k) {
      ulong const b_k = bm(k, S);
      ulong const S_b_k = S_n + b_k;
      ulong const S_b_nk = S_n + n * b_k;
      for (ulong j = 0; j < n; ++j) {
        if (j != b_k) {
          ulong const j_n = n * j;
          poly_t acc = field(0);
#pragma GCC unroll 4
#pragma GCC ivdep
          for (ulong l = 1; l <= bm_0_S; ++l) {
            auto const b_l = bm(l, S);
            acc += G0[j_n + b_l] * Xi[S_b_k + n * b_l];
            // Xi(b_k, b_l, S);
          }  // for l in range(0, n)
          rho_T[S_b_nk + j] = acc;
        }  // for k in range(0, n)
      }
    }  // for j in range(0, n)
  }    // for S in range(1, (1ul<<n))
  return;
}  // function CDet_rho_r

}  // namespace ffd::sigmadet_hubbard
