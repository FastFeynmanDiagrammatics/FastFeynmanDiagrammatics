namespace ffd::determinant {

template <typename matrix_t>

auto Gauss_CP_quad(matrix_t&& MM) {
  using ffd::nilpotent_polynomial::abs;
  using std::abs;
  using element_t =
      typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

  std::size_t const lin_size = ffd::core_math::sqrt_int(size(MM));
  assert(lin_size * lin_size == size(MM));

#if defined(__GNUC__) && !defined(__clang__)
  using quad_t = __float128;
#else
  using quad_t = long double;
#endif
  std::vector<quad_t> M(lin_size * lin_size);
  for (ulong j = 0; j < lin_size * lin_size; ++j) {
    M[j] = MM[j];
  }  // for j in range(0, lin_size*lin_size)

  switch (lin_size) {
    case 0ul:
      return element_t(1);
    case 1ul:
      return element_t(M[0]);
    case 2ul:
      return element_t(M[0] * M[3] - M[1] * M[2]);
    case 3ul:
      return element_t(M[0] * (M[4] * M[8] - M[5] * M[7]) -
                       M[1] * (M[3] * M[8] - M[5] * M[6]) +
                       M[2] * (M[3] * M[7] - M[4] * M[6]));
    default:
      quad_t det = 1.;

      for (std::size_t j = 0; j < lin_size - 1; ++j) {
        auto pivot_abs = abs<double>(M[j + lin_size * j]);
        std::size_t pivot_k = j;
        for (std::size_t k = j + 1; k < lin_size; ++k) {
          if (auto abs_temp = abs<double>(M[j + k * lin_size]);
              abs_temp > pivot_abs) {
            pivot_abs = abs_temp;
            pivot_k = k;
          }
        }
        bool swap_rows = true;
        for (std::size_t k = j + 1; k < lin_size; ++k) {
          if (auto abs_temp = abs<double>(M[k + j * lin_size]);
              abs_temp > pivot_abs) {
            pivot_abs = abs_temp;
            pivot_k = k;
            swap_rows = false;
          }
        }
        if (pivot_abs < 10000 * std::numeric_limits<quad_t>::min()) {
          return element_t(0.);
        }

        if (pivot_k != j) {
          if (swap_rows) {
            // swap rows
            for (std::size_t j2 = j; j2 < lin_size; ++j2) {
              std::swap(M[j2 + j * lin_size], M[j2 + pivot_k * lin_size]);
            }
            det = -det;
          } else {
            // swap columns
            for (std::size_t j2 = j; j2 < lin_size; ++j2) {
              std::swap(M[j + j2 * lin_size], M[pivot_k + j2 * lin_size]);
            }
            det = -det;
          }
        }

        // GAUSSIAN ELIMINATION
        det = M[j + j * lin_size] * det;
        auto const one_over_pivot = 1. / M[j + j * lin_size];
        for (std::size_t j2 = j + 1; j2 < lin_size; ++j2) {
          M[j2 + j * lin_size] = M[j2 + j * lin_size]
                                 // / M[j+j*lin_size];
                                 * one_over_pivot;
        }
        for (std::size_t k = j + 1; k < lin_size; ++k) {
          for (std::size_t j2 = j + 1; j2 < lin_size; ++j2) {
            M[j2 + k * lin_size] -= M[j + k * lin_size] * M[j2 + j * lin_size];
          }
        }
      }

      return element_t(det * M[lin_size - 1 + lin_size * (lin_size - 1)]);
  }
}

}  // namespace ffd::determinant
