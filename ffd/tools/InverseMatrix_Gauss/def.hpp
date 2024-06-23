namespace ffd::inverse_matrix_gauss {

template <typename matrix_t>

auto Inverse_and_Determinant(matrix_t&& M) {
  using ffd::nilpotent_polynomial::abs;
  using std::abs;
  using size_t = std::size_t;
  using element_t = typename std::decay_t<decltype(M[0])>;

  size_t const n = ffd::core_math::sqrt_int(size(M));
  assert(n * n == size(M));

  auto I = M;
  for (size_t j = 0; j < n; ++j) {
    for (size_t k = 0; k < n; ++k) {
      if (j != k) {
        I[k + j * n] = element_t{0.};
      } else {
        I[k + j * n] = element_t{1.};
      }
    }
  }

  element_t det{1.};

  switch (n) {
    case 0ul:
      break;
    case 1ul:
      det = M[0];
      assert((!assert_special_case || !NotZero(det)));
      I[0] = NotZero(det) ? element_t{1.} / det : element_t{0};
      break;
    case 2ul:
      det = M[0] * M[3] - M[1] * M[2];
      assert((!assert_special_case || !NotZero(det)));
      I[0] = NotZero(det) ? M[3] / det : element_t{0.};
      I[1] = NotZero(det) ? -M[1] / det : element_t{0.};
      I[2] = NotZero(det) ? -M[2] / det : element_t{0.};
      I[3] = NotZero(det) ? M[0] / det : element_t{0.};
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
          return std::make_pair(I, element_t(0.));
        }

        // SWAP ROWS
        if (pivot_k != j) {
          for (size_t j2 = j; j2 < n; ++j2) {
            std::swap(M[j2 + j * n], M[j2 + pivot_k * n]);
          }
          for (size_t j2 = 0; j2 < n; ++j2) {
            std::swap(I[j2 + j * n], I[j2 + pivot_k * n]);
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
          I[j2 + j * n] = I[j2 + j * n] * one_over_pivot;
        }
        for (size_t k = 0; k < j; ++k) {
          auto const M_kj = M[k * n + j];
          for (size_t j2 = 0; j2 < n; ++j2) {
            M[k * n + j2] -= M_kj * M[j * n + j2];
          }
          for (size_t j2 = 0; j2 < n; ++j2) {
            I[k * n + j2] -= M_kj * I[j * n + j2];
          }
        }
        for (size_t k = j + 1; k < n; ++k) {
          auto const M_kj = M[k * n + j];
          for (size_t j2 = 0; j2 < n; ++j2) {
            M[j2 + k * n] -= M_kj * M[j2 + j * n];
          }
          for (size_t j2 = 0; j2 < n; ++j2) {
            I[j2 + k * n] -= M_kj * I[j2 + j * n];
          }
        }
      }
  }

  return std::make_pair(I, det);
}

}  // namespace ffd::inverse_matrix_gauss
