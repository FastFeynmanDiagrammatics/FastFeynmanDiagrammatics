namespace ffd::canonical_sigma {

template <std::size_t n, class array_t>
auto ranked_zeta_transform(array_t const& f) {
  using el_t = std::decay_t<decltype(f[0])>;
  ffd::vector2d<el_t> ret((1ul << n), n + 1, 0.);
  ffd::vector2d<el_t> rzt((1ul << n), n + 1, 0.);
  for (ulong k = 0; k < n + 1; ++k) {
    for (ulong S = 0; S < (1ul << n); ++S) {
      rzt(S, 0) = __builtin_popcount(S) == k ? f[S] : 0;
    }  // for S in range(0, (1ul<<n))
    for (ulong j = 0; j < n; ++j) {
      for (ulong S = 0; S < (1ul << n); ++S) {
        if (((S >> j) & 1) == 0) {
          rzt(S, j + 1) = rzt(S, j);
        } else {
          rzt(S, j + 1) = rzt(S, j) + rzt(S - (1ul << j), j);
        }
      }  // for S in range(0, (1ul<<n))
    }    // for j in range(0, n+1)
    for (ulong S = 0; S < (1ul << n); ++S) {
      ret(S, k) = rzt(S, n);
    }  // for S in range(0, (1ul<<n))
  }    // for k in range(0, n+1)
  return ret;
}

// template <std::size_t n, class array_t>
// void ranked_zeta_transform_index_inv_order_r(
//     array_t& f,
//     std::array<std::size_t, 2> const k_array,
//     std::size_t const index = 0) {
//   using el_t = std::decay_t<decltype(f[0])>;
//   std::array<el_t, (1ul << n)> rzt;
//   for (ulong k = k_array[0]; k < k_array[1]; ++k) {
//     for (ulong S = 0; S < (1ul << n); ++S) {
//       rzt[S] = __builtin_popcount(S) == k ? f[index + S] : 0;
//     }  // for S in range(0, (1ul<<n))
//     for (ulong j = 0; j < n; ++j) {
//       for (ulong V = 0; V < (1ul << (n - j - 1)); ++V) {
//         auto const V_j = (V << (j + 1));
//         auto const V_j_p = V_j + (1ul << j);
//         for (ulong S = 0; S < (1ul << j); ++S) {
//           rzt[V_j_p + S] += rzt[V_j + S];
//         }  // for S
//       }    // for V
//     }      // for j
//     ulong const k_index = index + (k * (1ul << n);
//     for (ulong S = 0; S < (1ul << n); ++S) {
//       f[k_index + S] = rzt[S];
//     }  // for S in range(0, (1ul<<n))
//   }    // for k in range(0, n+1)
//   return;
// }

}  // namespace ffd::canonical_sigma
