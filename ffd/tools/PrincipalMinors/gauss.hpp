namespace ffd::principal_minors {

template <int order, int step = 1, class matrix_t = int>
auto gauss(matrix_t const& matrix) {
  using namespace ffd::nilpotent_polynomial;
  using element_t = typename std::decay_t<decltype(matrix[0])>;

  assert(order * order * step * step == (int)size(matrix));

  std::array<element_t, (1 << order)> minors;
  minors[0] = 1.;
  std::array<element_t, order * order * step * step> minor_matrix;
  std::array<ulong, order> bin_dig;
  for (ulong S = 1; S < (1ul << order); ++S) {
    ulong const card_S = __builtin_popcount(S);
    //  if (card_S > order / 2)
    // continue;
    ffd::set_theory::BinaryDigitsOf_r(S, bin_dig);
    if constexpr (step == 1) {
      for (ulong j = 0; j < card_S; ++j) {
        for (ulong k = 0; k < card_S; ++k) {
          minor_matrix[j * card_S + k] =
              matrix[bin_dig[j] * order + bin_dig[k]];
        }  // for k in range(0, card_S)
      }    // for j in range(0, card_S)
    } else if constexpr (step == 2) {
      for (ulong j = 0; j < card_S; ++j) {
        auto j_ind = 2 * j * 2 * card_S;
        auto bj_ind = 2 * bin_dig[j] * 2 * order;
        for (ulong k = 0; k < card_S; ++k) {
          minor_matrix[j_ind + 2 * k] = matrix[bj_ind + 2 * bin_dig[k]];
          minor_matrix[j_ind + 2 * k + 1] = matrix[bj_ind + 2 * bin_dig[k] + 1];
          minor_matrix[j_ind + 2 * card_S + 2 * k] =
              matrix[bj_ind + 2 * order + 2 * bin_dig[k]];
          minor_matrix[j_ind + 2 * card_S + 2 * k + 1] =
              matrix[bj_ind + 2 * order + 2 * bin_dig[k] + 1];
        }  // for k in range(0, card_S)
      }    // for j in range(0, card_S)
    }
    minors[S] = ffd::determinant::gauss(minor_matrix, step * card_S);
  }
  return minors;
}

}  // namespace ffd::principal_minors
