namespace ffd::principal_minors {

template <class high_prec_t = long double, typename matrix_t = int>
auto PrincipalMinors2n_nopivot_high(matrix_t&& matrix) {
  using element_t =
      typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

  using namespace ffd::nilpotent_polynomial;
  // using ffd::vector_range::Range;
  using ffd::math_tools::print_matrix;
  using std::abs;
  using std::cerr;
  using std::cout;
  using std::endl;

  const BinaryInt order = ffd::core_math::sqrt_int(size(matrix)) / 2;
  assert(4 * order * order == (int)size(matrix));
  const BinaryInt card = (1 << order);

  const BinaryInt max_list_size = 4 * 9 * (1 << ((order - 3) * (order > 2)));

  std::vector<high_prec_t> matrix_list_old(max_list_size);
  std::vector<high_prec_t> matrix_list_new(max_list_size);

  for (BinaryInt element = 0; element < size(matrix); ++element) {
    matrix_list_old[element] = matrix[element];
  }

  std::vector<high_prec_t> minors(card);
  minors[0] = high_prec_t(1);
  minors[1] += matrix_list_old[0];
  BinaryInt mnr = 1;

  std::vector<high_prec_t> inverse(4);
  std::vector<high_prec_t> inverse_1(4 * order);

  // BREADTH FIRST TRAVERSAL
  for (BinaryInt level = 0; level < order - 1; ++level) {
    const BinaryInt old_lin_size = 2 * (order - level);
    const BinaryInt old_mat_size = old_lin_size * old_lin_size;
    const BinaryInt new_lin_size = old_lin_size - 2;
    const BinaryInt new_mat_size = new_lin_size * new_lin_size;
    const BinaryInt two_level = (1 << level);
    const BinaryInt shift = two_level * new_mat_size;
    for (BinaryInt set = 0; set < two_level; ++set) {
      BinaryInt pos1 = set * new_mat_size;
      BinaryInt pos2 = set * old_mat_size;
      BinaryInt pos3 = pos2 + 2 * old_lin_size;
      BinaryInt pos3_ = pos3 + 1;
      BinaryInt pos4 = pos3 + 2;

      inverse[0] = matrix_list_old[pos2 + old_lin_size + 1];
      inverse[1] = -matrix_list_old[pos2 + 1];
      inverse[2] = -matrix_list_old[pos2 + old_lin_size];
      inverse[3] = matrix_list_old[pos2];

      const high_prec_t determinant =
          inverse[0] * inverse[3] - inverse[1] * inverse[2];
      inverse[0] /= determinant;
      inverse[1] /= determinant;
      inverse[2] /= determinant;
      inverse[3] /= determinant;

      minors[mnr++] = determinant;

#pragma unroll
      for (BinaryInt j = 0; j < new_lin_size; ++j) {
        const high_prec_t m1 = matrix_list_old[pos2 + 2 + j];
        const high_prec_t m2 = matrix_list_old[pos2 + 2 + j + old_lin_size];
        inverse_1[j] = m1 * inverse[0] + m2 * inverse[1];
        inverse_1[j + new_lin_size] = m1 * inverse[2] + m2 * inverse[3];
      }

#pragma unroll
      for (BinaryInt i = 0; i < new_lin_size; ++i) {
#pragma unroll
        for (BinaryInt j = 0; j < new_lin_size; ++j) {
          const high_prec_t el = matrix_list_old[pos4];
          matrix_list_new[pos1] = el;
          matrix_list_new[pos1 + shift] =
              el - (matrix_list_old[pos3] * inverse_1[j] +
                    matrix_list_old[pos3_] * inverse_1[j + new_lin_size]);
          ++pos1;
          ++pos4;
        }
        pos4 += 2;
        pos3 += old_lin_size;
        pos3_ += old_lin_size;
      }
    }
    // print_matrix<2>(matrix_list_new, new_lin_size*(1<<(level+1)),
    // new_lin_size);
    matrix_list_old.swap(matrix_list_new);
  }

  for (BinaryInt set = 0; set < (1 << (order - 1)); ++set) {
    const BinaryInt pos2 = set * 4;
    minors[mnr++] = matrix_list_old[pos2] * matrix_list_old[pos2 + 3] -
                    matrix_list_old[pos2 + 1] * matrix_list_old[pos2 + 2];
    // if (abs(determinant) < epsilon)
    // return DeterminantPrincipalMinors(matrix, 2);
  }

  // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
  BinaryInt set_K = 1;
  for (BinaryInt level = 0; level < order; ++level) {
    const BinaryInt set_2K = 2 * set_K;
    for (BinaryInt set = set_K + 1; set < set_2K; ++set) {
      minors[set] *= minors[set - set_K];
    }
    set_K = set_2K;
  }

  std::vector<element_t> minors_ret(card);
  for (ulong j = 0; j < minors_ret.size(); ++j) {
    minors_ret[j] = minors[j];
  }  // for j in range(0, minors_ret.size())

  return minors_ret;
}

}  // namespace ffd::principal_minors
