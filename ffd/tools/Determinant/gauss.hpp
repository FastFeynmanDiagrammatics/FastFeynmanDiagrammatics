namespace ffd::determinant {

template <typename matrix_t>

auto gauss(matrix_t& M, ulong const n) {
  using ffd::nilpotent_polynomial::abs;
  using std::abs;
  using element_t = typename std::decay_t<decltype(M[0])>;

  switch (n) {
    case 0ul:
      return element_t(1.);
    case 1ul:
      return M[0];
    case 2ul:
      return M[0] * M[3] - M[1] * M[2];
    case 3ul:
      return M[0] * (M[4] * M[8] - M[5] * M[7]) -
             M[1] * (M[3] * M[8] - M[5] * M[6]) +
             M[2] * (M[3] * M[7] - M[4] * M[6]);
    case 4ul:
      return +(M[0] * M[5] - M[4] * M[1]) * (M[10] * M[15] - M[11] * M[14]) -
             (M[0] * M[6] - M[4] * M[2]) * (M[9] * M[15] - M[11] * M[13]) +
             (M[0] * M[7] - M[4] * M[3]) * (M[9] * M[14] - M[10] * M[13]) +
             (M[1] * M[6] - M[5] * M[2]) * (M[8] * M[15] - M[11] * M[12]) -
             (M[1] * M[7] - M[5] * M[3]) * (M[8] * M[14] - M[10] * M[12]) +
             (M[2] * M[7] - M[6] * M[3]) * (M[8] * M[13] - M[9] * M[12]);
    default:
      element_t det = 1.;

      for (std::size_t j = 0; j < n - 4; ++j) {
        auto pivot_abs = abs(M[j + j * n]);
        std::size_t pivot_k = j;
        for (std::size_t k = j + 1; k < n; ++k) {
          if (auto const abs_temp = abs(M[j + k * n]); abs_temp > pivot_abs) {
            pivot_abs = abs_temp;
            pivot_k = k;
          }
        }
        if (pivot_abs < 10000 * std::numeric_limits<Real>::min()) {
          return element_t(0.);
        }

        // SWAP ROWS
        if (pivot_k != j) {
          for (std::size_t j2 = j; j2 < n; ++j2) {
            std::swap(M[j2 + j * n], M[j2 + pivot_k * n]);
          }
          det = -det;
        }

        // GAUSSIAN ELIMINATION
        det = M[j + j * n] * det;
        auto const one_over_pivot = 1. / M[j + j * n];
        for (std::size_t j2 = j + 1; j2 < n; ++j2) {
          M[j2 + j * n] = M[j2 + j * n]
                          // / M[j+j*n];
                          * one_over_pivot;
        }
        //        for (std::size_t k_n = (j + 1) * n; k_n < n * n; k_n += n) {
        for (std::size_t k = j + 1; k < n; k++) {
          for (std::size_t j2 = j + 1; j2 < n; ++j2) {
            M[j2 + k * n] -= M[j + k * n] * M[j2 + j * n];
          }
        }
      }

      return det * (+(M[n - 4 + (n - 4) * n] * M[n - 3 + n * (n - 3)] -
                      M[n - 4 + n * (n - 3)] * M[n - 3 + n * (n - 4)]) *
                        (M[n - 2 + n * (n - 2)] * M[n - 1 + n * (n - 1)] -
                         M[n - 1 + n * (n - 2)] * M[n - 2 + n * (n - 1)]) -
                    (M[n - 4 + n * (n - 4)] * M[n - 2 + n * (n - 3)] -
                     M[n - 4 + n * (n - 3)] * M[n - 2 + n * (n - 4)]) *
                        (M[n - 3 + n * (n - 2)] * M[n - 1 + n * (n - 1)] -
                         M[n - 1 + n * (n - 2)] * M[n - 3 + n * (n - 1)]) +
                    (M[n - 4 + n * (n - 4)] * M[n - 1 + n * (n - 3)] -
                     M[n - 4 + n * (n - 3)] * M[n - 1 + n * (n - 4)]) *
                        (M[n - 3 + n * (n - 2)] * M[n - 2 + n * (n - 1)] -
                         M[n - 2 + n * (n - 2)] * M[n - 3 + n * (n - 1)]) +
                    (M[n - 3 + n * (n - 4)] * M[n - 2 + n * (n - 3)] -
                     M[n - 3 + n * (n - 3)] * M[n - 2 + n * (n - 4)]) *
                        (M[n - 4 + n * (n - 2)] * M[n - 1 + n * (n - 1)] -
                         M[n - 1 + n * (n - 2)] * M[n - 4 + n * (n - 1)]) -
                    (M[n - 3 + n * (n - 4)] * M[n - 1 + n * (n - 3)] -
                     M[n - 3 + n * (n - 3)] * M[n - 1 + n * (n - 4)]) *
                        (M[n - 4 + n * (n - 2)] * M[n - 2 + n * (n - 1)] -
                         M[n - 2 + n * (n - 2)] * M[n - 4 + n * (n - 1)]) +
                    (M[n - 2 + n * (n - 4)] * M[n - 1 + n * (n - 3)] -
                     M[n - 2 + n * (n - 3)] * M[n - 1 + n * (n - 4)]) *
                        (M[n - 4 + n * (n - 2)] * M[n - 3 + n * (n - 1)] -
                         M[n - 3 + n * (n - 2)] * M[n - 4 + n * (n - 1)]));
  }
}

}  // namespace ffd::determinant
