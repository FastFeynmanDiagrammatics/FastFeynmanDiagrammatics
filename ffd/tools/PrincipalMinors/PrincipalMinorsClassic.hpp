namespace ffd::principal_minors {

template <typename matrix_t>
std::vector<typename std::decay<decltype(std::declval<matrix_t>()[0])>::type>
PrincipalMinorsBFS(matrix_t&& matrix) {
  using element_t =
      typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
  using namespace ffd::nilpotent_polynomial;
  // using ffd::vector_range::Range;
  using std::abs;
  // using std::cerr;
  // using std::cout;
  // using std::endl;
  // using ffd::math_tools::print_matrix;

  const BinaryInt order = ffd::core_math::sqrt_int(size(matrix));
  assert(order * order == (int)size(matrix));
  const BinaryInt card = (1 << order);

  element_t pivot_ = 0.0;
  for (BinaryInt element = 0; element < size(matrix); ++element) {
    pivot_ = std::max(abs(matrix[element]), pivot_);
  }
  pivot_ /= (Real)size(matrix);
  const element_t pivot = pivot_;

  const BinaryInt max_list_size = 9 * (1 << ((order - 3) * (order > 2)));

  std::vector<element_t> matrix_list_old(max_list_size);
  std::vector<element_t> matrix_list_new(max_list_size);
  // std::vector<element_t> matrix_list_old;
  // std::vector<element_t> matrix_list_new;
  // matrix_list_old.reserve(max_list_size);
  // matrix_list_new.reserve(max_list_size);

  for (BinaryInt element = 0; element < size(matrix); ++element) {
    matrix_list_old[element] = matrix[element];
  }

  std::vector<element_t> minors(card, pivot);
  minors[0] = 1.;
  minors[1] += matrix_list_old[0];
  BinaryInt mnr = 2;

  // print_matrix<2>(matrix_list_old, order, order);

  // BREADTH FIRST TRAVERSAL
  for (BinaryInt level = 0; level < order - 1; ++level) {
    const BinaryInt old_lin_size = order - level;
    const BinaryInt old_mat_size = old_lin_size * old_lin_size;
    const BinaryInt new_lin_size = old_lin_size - 1;
    const BinaryInt new_mat_size = new_lin_size * new_lin_size;
    const BinaryInt two_level = (1 << level);
    for (BinaryInt set = 0; set < two_level; ++set) {
      // std::cout << "set " << set << std::endl;
      BinaryInt pos1 = set * new_mat_size;
      BinaryInt pos2 = set * old_mat_size;
      BinaryInt pos3 = pos2 + old_lin_size;
      BinaryInt pos4 = pos3 + 1;
      BinaryInt pos5 = pos1 + two_level * new_mat_size;
      const element_t inverse = 1.0 / (matrix_list_old[pos2] + pivot);
#pragma unroll
      for (BinaryInt i = 0; i < new_lin_size; ++i) {
        BinaryInt pos6 = pos2 + 1;
        const element_t inverse_1 = matrix_list_old[pos3] * inverse;
#pragma unroll
        for (BinaryInt j = 0; j < new_lin_size; ++j) {
          const element_t el = matrix_list_old[pos4];
          matrix_list_new[pos1] = el;
          matrix_list_new[pos5] = el - matrix_list_old[pos6] * inverse_1;
          ++pos1;
          ++pos4;
          ++pos5;
          ++pos6;
        }
        ++pos4;
        pos3 += old_lin_size;
      }
    }
    for (BinaryInt set = 0; set < two_level * 2 * new_mat_size;
         set += new_mat_size) {
      // std::cout <<
      minors[mnr++] += matrix_list_new[set];
    }
    // print_matrix<2>(matrix_list_new, new_lin_size*(1<<(level+1)),
    // new_lin_size);
    matrix_list_old.swap(matrix_list_new);
  }

  // std::cout << "after big loop:" << std::endl;
  // std::cout << std::setprecision(10);
  // for (uint m=0; m<(1<<order); ++m){
  //   std::cout << m << " " << minors[m] << std::endl;
  // }

  // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
  BinaryInt set_K = 1;
  for (BinaryInt level = 0; level < order; ++level) {
    const BinaryInt set_2K = 2 * set_K;
    for (BinaryInt set = set_K + 1; set < set_2K; ++set) {
      minors[set] *= minors[set - set_K];
    }
    set_K = set_2K;
  }

  // std::cout << "after multiply loop:" << std::endl;
  // for (uint m=0; m<(1<<order); ++m){
  //   std::cout << m << " " << minors[m] << std::endl;
  // }

  // PIVOT CORRECTIONS
  set_K = card;
  for (BinaryInt level = order; level > 0; --level) {
    const BinaryInt set_K_2 = (set_K >> 1);
    for (BinaryInt set = set_K - 1; set >= set_K_2; --set) {
      for (BinaryInt set_ = set; set_ < card; set_ += set_K) {
        minors[set_] -= pivot * minors[set_ - set_K_2];
      }
    }
    set_K = set_K_2;
  }

  // std::cout << "after pivot loop:" << std::endl;
  // for (uint m=0; m<(1<<order); ++m){
  //   std::cout << m << " " << minors[m] << std::endl;
  // }

  return minors;
}

template <typename matrix_t>
std::vector<typename std::decay<decltype(std::declval<matrix_t>()[0])>::type>
PrincipalMinorsBFS2n(matrix_t&& matrix, double epsilon = 1.e-5) {
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

  std::vector<element_t> matrix_list_old(max_list_size);
  std::vector<element_t> matrix_list_new(max_list_size);

  for (BinaryInt element = 0; element < size(matrix); ++element) {
    matrix_list_old[element] = matrix[element];
  }

  std::vector<element_t> minors(card);
  minors[0] = 1.;
  minors[1] += matrix_list_old[0];
  BinaryInt mnr = 1;

  std::vector<element_t> inverse(4);
  std::vector<element_t> inverse_1(4 * order);

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

      const element_t determinant =
          inverse[0] * inverse[3] - inverse[1] * inverse[2];
      if (abs(determinant) < epsilon) {
        std::cout << "going to dets" << std::endl;
        return DeterminantPrincipalMinors(matrix, 2);
      }
      inverse[0] /= determinant;
      inverse[1] /= determinant;
      inverse[2] /= determinant;
      inverse[3] /= determinant;

      minors[mnr++] = determinant;

#pragma unroll
      for (BinaryInt j = 0; j < new_lin_size; ++j) {
        const element_t m1 = matrix_list_old[pos2 + 2 + j];
        const element_t m2 = matrix_list_old[pos2 + 2 + j + old_lin_size];
        inverse_1[j] = m1 * inverse[0] + m2 * inverse[1];
        inverse_1[j + new_lin_size] = m1 * inverse[2] + m2 * inverse[3];
      }

#pragma unroll
      for (BinaryInt i = 0; i < new_lin_size; ++i) {
#pragma unroll
        for (BinaryInt j = 0; j < new_lin_size; ++j) {
          const element_t el = matrix_list_old[pos4];
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

  return minors;
}

template <typename matrix_t>
std::vector<typename std::decay<decltype(std::declval<matrix_t>()[0])>::type>
PrincipalMinorsPivot2n(matrix_t&& matrix) {
  using element_t =
      typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
  using std::abs;

  const BinaryInt order = ffd::core_math::sqrt_int(size(matrix)) / 2;
  assert(4 * order * order == (int)size(matrix));
  const BinaryInt card = (1 << order);

  const BinaryInt max_list_size = 4 * 9 * (1 << ((order - 3) * (order > 2)));

  std::vector<element_t> matrix_list_old(max_list_size);
  std::vector<element_t> matrix_list_new(max_list_size);

  for (BinaryInt element = 0; element < size(matrix); ++element) {
    matrix_list_old[element] = matrix[element];
  }

  std::vector<element_t> inverse(4);
  std::vector<element_t> inverse_1(4 * order);
  std::vector<element_t> minors(card, 0.);

  // compute zeroth minor
  minors[0] = 1.;
  minors[1] += matrix_list_old[0] * matrix_list_old[1 + 2 * order] -
               matrix_list_old[1] * matrix_list_old[2 * order];
  BinaryInt mnr = 2;

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

      const element_t determinant =
          inverse[0] * inverse[3] - inverse[1] * inverse[2];

      inverse[0] /= determinant;
      inverse[1] /= determinant;
      inverse[2] /= determinant;
      inverse[3] /= determinant;

#pragma unroll
      for (BinaryInt j = 0; j < new_lin_size; ++j) {
        const element_t m1 = matrix_list_old[pos2 + 2 + j];
        const element_t m2 = matrix_list_old[pos2 + 2 + j + old_lin_size];
        inverse_1[j] = m1 * inverse[0] + m2 * inverse[1];
        inverse_1[j + new_lin_size] = m1 * inverse[2] + m2 * inverse[3];
      }

#pragma unroll
      for (BinaryInt i = 0; i < new_lin_size; ++i) {
#pragma unroll
        for (BinaryInt j = 0; j < new_lin_size; ++j) {
          const element_t el = matrix_list_old[pos4];
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

    for (BinaryInt set = 0; set < two_level * 2 * new_mat_size;
         set += new_mat_size) {
      minors[mnr++] +=
          matrix_list_new[set] * matrix_list_new[set + 1 + new_lin_size] -
          matrix_list_new[set + 1] * matrix_list_new[set + new_lin_size];
    }
    matrix_list_old.swap(matrix_list_new);
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

  return minors;
}

}  // namespace ffd::principal_minors