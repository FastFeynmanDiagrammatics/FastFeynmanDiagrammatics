namespace ffd::canonical_sigma {

// template <std::size_t n, std::size_t poly_s, class field>
// void Xi_r(array2d<Real, n, n> const& __restrict__ G0,
//           vector3d<field>& __restrict__ Xi,
//           array2d<std::uint8_t, n + 1, (1ul << n)> const& __restrict__ bm) {
//   //    t_3.ini();

//   vector4d<Real> xi((1ul << n), poly_s + 1, n, n, Real(0));
//   vector2d<Real> z((1ul << n), poly_s + 1, Real(0));
//   vector1d<field> Z(1ul << n);
//   z(0, 0) = 1;
//   for (BinaryInt S = 1; S < (1ul << n); ++S) {
//     ulong const card_S = bm(0, S);
//     bool const is_even_S = ((card_S & 1) == 0);

//     vector2d<Real> g0(card_S, card_S);
//     //  t_2.ini();
//     for (ulong j = 0; j < card_S; ++j) {
//       auto const bm_j1_S = bm(j + 1, S);
//       for (ulong k = 0; k < card_S; ++k) {
//         g0(k, j) = G0(bm(k + 1, S), bm_j1_S);
//       }
//     }
//     // t_2.fin();
//     //      t_temp.ini();
//     auto const [I, det] =
//         ffd::inverse_matrix_gauss::Inverse_and_Determinant(std::move(g0));
//     //      t_temp.fin();
//     z(S, 0) = det;

//     for (ulong j = 0; j < card_S; ++j) {
//       auto const bm_j1_S = bm(j + 1, S);
//       for (ulong k = 0; k < card_S; ++k) {
//         auto const bm_k1_S = bm(k + 1, S);
//         if (bm_k1_S != bm_j1_S)
//           xi(S - (1ul << bm_j1_S) - (1ul << bm_k1_S), 0, bm(k + 1, S),
//              bm_j1_S) = I[j * card_S + k] * det;
//         else
//           xi(S - (1ul << bm_j1_S), 0, bm(k + 1, S), bm_j1_S) =
//               I[j * card_S + k] * det;
//       }
//     }
//   }

//   ranked_zeta_transform_index_r<n>(z, poly_s);
//   ulong const size_rzt = (1ul << n) * (poly_s + 1);
//   for (ulong j = 0; j < n; ++j) {
//     ulong const j_index = size_rzt * n * j;
//     for (ulong k = 0; k < n; ++k) {
//       ranked_zeta_transform_index_r<n>(xi, poly_s, j_index + size_rzt * k);
//     }  // for k in range(0, n)
//   }    // for j in range(0, n)// compute connected part
//   //    t_1.ini();

//   for (ulong S = 0; S < (1ul << n); ++S) {
//     auto const card_S = bm(0, S);
//     for (ulong p = 0; p < poly_s + 1; ++p) {
//       Z[S][p] = z(S, p);
//     }
//     if (card_S >= 2) {
//       for (ulong j = 0; j < card_S; ++j) {
//         auto const bm_j = bm(j + 1, S);
//         auto const S_j = S - (1ul << bm_j);
//         for (ulong k = 0; k < card_S; ++k) {
//           auto const bm_k = bm(k + 1, S);
//           if (j != k) {
//             auto const S_j_k = S_j - (1ul << bm_k);
//             for (ulong p = 0; p < poly_s + 1; ++p) {
//               Xi(bm_k, bm_j, S)[p] = xi(S_j_k, p, bm_k, bm_j);
//             }  // for p in range(0, poly_s+1)
//           } else {
//             for (ulong p = 0; p < poly_s + 1; ++p) {
//               Xi(bm_k, bm_j, S)[p] = xi(S_j, p, bm_k, bm_j);
//             }  // for p in range(0, poly_s+1)
//           }
//         }  // for k in range(0, n)
//       }    // for j in range(0, n)
//     }      // if card_S >= 2
//   }        // for S in range(0, (1ul<<n))

//   for (ulong S = 0; S < (1ul << n); ++S) {
//     auto const card_S = bm(0, S);
//     for (ulong j = 0; j < card_S; ++j) {
//       auto const bm_j = bm(j + 1, S);
//       for (ulong k = 0; k < card_S; ++k) {
//         auto const bm_k = bm(k + 1, S);
//         Xi(bm_k, bm_j, S) *= Z[S];
//         if ((card_S & 1) == 0) {
//           Xi(bm_k, bm_j, S) = -Xi(bm_k, bm_j, S);
//         }  // if
//       }    // for k
//     }      // for j
//     Z[S] *= Z[S];
//     if ((card_S & 1) != 0) {
//       Z[S] = -Z[S];
//     }
//   }  // for S

//   ulong constexpr n2 = n * n;
//   for (ulong V = 1; V < (1ul << n); ++V) {
//     auto const V_n = n2 * V;
//     for (ulong S = ((V - 1) & V); S != 0; S = ((S - 1) & V)) {
//       auto const S_n = n2 * S;
//       auto const Z_V_S = Z[V - S];
//       ulong const bm_0_S = bm(0, S);
//       for (ulong j = 0; j < bm_0_S; ++j) {
//         ulong const bmj_n = bm(j + 1, S) * n;
//         ulong const V_j_n = V_n + bmj_n;
//         ulong const S_j_n = S_n + bmj_n;
// #pragma GCC unroll 4
// #pragma GCC ivdep
//         for (ulong k = 0; k < bm_0_S; ++k) {
//           auto const bm_k1_S = bm(k + 1, S);
//           Xi[V_j_n + bm_k1_S] -= Xi[S_j_n + bm_k1_S] * Z_V_S;
//         }
//       }  // for j in range(0, bm(0, S))
//     }
//   }
//   // t_1.fin();
//   // t_3.fin();
//   return;
// }  // namespace ffd::canonical_sigma

}  // namespace ffd::canonical_sigma
