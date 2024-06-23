namespace ffd::principal_minors {
  
  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  PrincipalMinors2n_nopivot(matrix_t&& matrix){
    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
    using std::abs;
    const BinaryInt order = ffd::core_math::sqrt_int(size(matrix))/2;
    assert (4*order*order == (int)size(matrix));
    const BinaryInt card = (1 << order);
    const BinaryInt max_list_size=4*9*(1<<((order-3)*(order>2)));
    std::vector<element_t> matrix_list_old(max_list_size);
    std::vector<element_t> matrix_list_new(max_list_size);
    for (BinaryInt element = 0; element < size(matrix); ++element){
      matrix_list_old[element] = matrix[element];
    }
    std::vector<element_t> inverse(4);
    std::vector<element_t> inverse_1(4*order);
    std::vector<element_t> minors(card,0.);
    // compute zeroth minor
    minors[0] = 1.;
    minors[1] += matrix_list_old[0]*matrix_list_old[1+2*order]
                -matrix_list_old[1]*matrix_list_old[2*order];
    BinaryInt mnr=2;
    // BREADTH FIRST TRAVERSAL
    for (BinaryInt level=0; level<order-1; ++level){
      const BinaryInt old_lin_size = 2*(order - level);
      const BinaryInt old_mat_size = old_lin_size * old_lin_size;
      const BinaryInt new_lin_size = old_lin_size - 2;
      const BinaryInt new_mat_size = new_lin_size * new_lin_size;
      const BinaryInt two_level    = (1<<level);
      const BinaryInt shift        = two_level * new_mat_size;
      for (BinaryInt set=0; set<two_level; ++set){
        BinaryInt pos1 = set * new_mat_size;
        BinaryInt pos2 = set * old_mat_size;
        BinaryInt pos3 = pos2 + 2*old_lin_size;
        BinaryInt pos3_= pos3 + 1;
        BinaryInt pos4 = pos3 + 2;
        inverse[0b00] =  matrix_list_old[pos2+old_lin_size+1];
        inverse[0b01] = -matrix_list_old[pos2+1];
        inverse[0b10] = -matrix_list_old[pos2+old_lin_size];
        inverse[0b11] =  matrix_list_old[pos2];
        const element_t determinant_0 = inverse[0b00] - inverse[0b01]*(inverse[0b10]/inverse[0b11]);
        inverse[0b00] /= determinant_0; inverse[0b00] /= inverse[0b11];
        inverse[0b01] /= determinant_0; inverse[0b01] /= inverse[0b11];
        inverse[0b10] /= determinant_0; inverse[0b10] /= inverse[0b11];
        inverse[0b11] /= inverse[0b11]; inverse[0b11] /= determinant_0;
        #pragma unroll
        for (BinaryInt j=0; j<new_lin_size; ++j){
          const element_t m1 = matrix_list_old[pos2+2+j];
          const element_t m2 = matrix_list_old[pos2+2+j+old_lin_size];
          inverse_1[j]              = m1 * inverse[0b00] + m2 * inverse[0b01];
          inverse_1[j+new_lin_size] = m1 * inverse[0b10] + m2 * inverse[0b11];
        }
        #pragma unroll
        for (BinaryInt i = 0; i<new_lin_size; ++i){
          #pragma unroll
          for (BinaryInt j = 0; j<new_lin_size; ++j){
            const element_t el = matrix_list_old[pos4];
            matrix_list_new[pos1] = el;
            matrix_list_new[pos1+shift] = el -
               (matrix_list_old[pos3] * inverse_1[j]
              + matrix_list_old[pos3_] * inverse_1[j+new_lin_size]);
            ++pos1;
            ++pos4;
          }
          pos4+=2;
          pos3+=old_lin_size;
          pos3_+=old_lin_size;
        }
      }
      for (BinaryInt set=0; set<two_level*2*new_mat_size; set+=new_mat_size){
        minors[mnr++] += matrix_list_new[set]*matrix_list_new[set+1+new_lin_size]
                       - matrix_list_new[set+1]*matrix_list_new[set+new_lin_size];
      }
      matrix_list_old.swap(matrix_list_new);
    }
    // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
    BinaryInt set_K = 1;
    for (BinaryInt level = 0; level < order; ++level){
      const BinaryInt set_2K = 2*set_K;
      for (BinaryInt set = set_K+1; set < set_2K; ++set){
        minors[set] *= minors[set-set_K];
      }
      set_K=set_2K;
    }
    return minors;
  }


} // namespace
