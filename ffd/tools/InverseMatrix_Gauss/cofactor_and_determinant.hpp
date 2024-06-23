namespace ffd::inverse_matrix_gauss {

template <typename matrix_t>

auto Cofactor_and_Determinant(matrix_t&& M) {
  using ffd::nilpotent_polynomial::abs;
  using std::abs;
  using size_t = std::size_t;
  using element_t = typename std::decay_t<decltype(M[0])>;

  size_t const n = ffd::core_math::sqrt_int(size(M));
  assert(n * n == size(M));

  auto C = M;
  for (size_t j = 0; j < n; ++j) {
    for (size_t k = 0; k < n; ++k) {
      if (j != k) {
        C[k + j * n] = element_t{0.};
      } else {
        C[k + j * n] = element_t{1.};
      }
    }
  }

  element_t det{1.};

  switch (n) {
    case 0ul:
      break;
    case 1ul:
      det = M[0];
      break;
    case 2ul:
      det = M[0b00] * M[0b11] - M[0b01] * M[0b10];
      C[0b00] = M[0b11];
      C[0b01] = -M[0b01];
      C[0b10] = -M[0b10];
      C[0b11] = M[0b00];
      break;
    case 3ul:
      det = M[0] * (M[4] * M[8] - M[5] * M[7]) -
            M[1] * (M[3] * M[8] - M[5] * M[6]) +
            M[2] * (M[3] * M[7] - M[4] * M[6]);
      C[0] = (M[4] * M[8] - M[5] * M[7]);
      C[1] = -(M[1] * M[8] - M[2] * M[7]);
      C[2] = (M[1] * M[5] - M[2] * M[4]);
      C[3] = -(M[3] * M[8] - M[5] * M[6]);
      C[4] = (M[0] * M[8] - M[2] * M[6]);
      C[5] = -(M[0] * M[5] - M[2] * M[3]);
      C[6] = (M[3] * M[7] - M[4] * M[6]);
      C[7] = -(M[0] * M[7] - M[1] * M[6]);
      C[8] = (M[0] * M[4] - M[1] * M[3]);
      break;
    default:
      for (size_t j = 0; j < n; ++j) {
        auto pivot_abs = abs(M[j + n * j]);
        size_t pivot_k = j;
        for (size_t k = j + 1; k < n; ++k) {
          if (abs(M[j + k * n]) > pivot_abs) {
            pivot_abs = abs(M[j + k * n]);
            pivot_k = k;
          }
        }
        if (pivot_abs < 10000 * std::numeric_limits<Real>::min()) {
          return std::make_pair(C, element_t(0.));
        }

        // SWAP ROWS
        if (pivot_k != j) {
          for (size_t j2 = j; j2 < n; ++j2) {
            std::swap(M[j2 + j * n], M[j2 + pivot_k * n]);
          }
          for (size_t j2 = 0; j2 < n; ++j2) {
            std::swap(C[j2 + j * n], C[j2 + pivot_k * n]);
          }
          det = -det;
        }

        // GAUSSIAN ELIMINATION
        det = M[j + j * n] * det;
        auto const one_over_pivot = 1. / M[j + j * n];
        for (size_t j2 = 0; j2 < n; ++j2) {
          M[j * n + j2] = M[j * n + j2] * one_over_pivot;
        }
        for (size_t j2 = 0; j2 < n; ++j2) {
          C[j2 + j * n] = C[j2 + j * n] * one_over_pivot;
        }

        for (size_t k = 0; k < j; ++k) {
          auto const M_kj = M[k * n + j];
          for (size_t j2 = 0; j2 < n; ++j2) {
            M[k * n + j2] -= M_kj * M[j * n + j2];
          }
          for (size_t j2 = 0; j2 < n; ++j2) {
            C[k * n + j2] -= M_kj * C[j * n + j2];
          }
        }
        for (size_t k = j + 1; k < n; ++k) {
          auto const M_kj = M[k * n + j];
          for (size_t j2 = 0; j2 < n; ++j2) {
            M[j2 + k * n] -= M_kj * M[j2 + j * n];
          }
          for (size_t j2 = 0; j2 < n; ++j2) {
            C[j2 + k * n] -= M_kj * C[j2 + j * n];
          }
        }
      }

      for (ulong u = 0; u < n * n; ++u) {
        C[u] *= det;
      }  // for u in range(0, n*n)

  }  // switch n

  return std::make_pair(C, det);
}

}  // namespace ffd::inverse_matrix_gauss
