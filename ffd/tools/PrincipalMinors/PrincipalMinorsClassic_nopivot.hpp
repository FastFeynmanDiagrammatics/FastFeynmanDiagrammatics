namespace ffd::principal_minors{

  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  PrincipalMinorsBFS_nopivot(matrix_t&& matrix){

    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
    using namespace ffd::nilpotent_polynomial;
    // using ffd::vector_range::Range;
    using std::abs;
    // using std::cerr;
    // using std::cout;
    // using std::endl;
    // using ffd::math_tools::print_matrix;

    const BinaryInt order = ffd::core_math::sqrt_int(size(matrix));
    assert (order*order == (int)size(matrix));
    const BinaryInt card = (1 << order);

    const BinaryInt max_list_size=9*(1<<((order-3)*(order>2)));

    std::vector<element_t> matrix_list_old(max_list_size);
    std::vector<element_t> matrix_list_new(max_list_size);
    // std::vector<element_t> matrix_list_old;
    // std::vector<element_t> matrix_list_new;
    // matrix_list_old.reserve(max_list_size);
    // matrix_list_new.reserve(max_list_size);

    for (BinaryInt element = 0; element < size(matrix); ++element){
      matrix_list_old[element] = matrix[element];
    }

    std::vector<element_t> minors(card, 0.);
    minors[0] = 1.;
    minors[1] += matrix_list_old[0];
    BinaryInt mnr=2;

    // print_matrix<2>(matrix_list_old, order, order);

    // BREADTH FIRST TRAVERSAL
    for (BinaryInt level=0; level<order-1; ++level){
      const BinaryInt old_lin_size = order - level;
      const BinaryInt old_mat_size = old_lin_size * old_lin_size;
      const BinaryInt new_lin_size = old_lin_size - 1;
      const BinaryInt new_mat_size = new_lin_size * new_lin_size;
      const BinaryInt two_level    = (1<<level);
      for (BinaryInt set=0; set<two_level; ++set){
        BinaryInt pos1 = set * new_mat_size;
        BinaryInt pos2 = set * old_mat_size;
        BinaryInt pos3 = pos2 + old_lin_size;
        BinaryInt pos4 = pos3 + 1;
        BinaryInt pos5 = pos1 + two_level * new_mat_size;
        const element_t inverse = 1.0/matrix_list_old[pos2];
        #pragma unroll
        for (BinaryInt i = 0; i<new_lin_size; ++i){
          BinaryInt pos6 = pos2 + 1;
          const element_t inverse_1 = matrix_list_old[pos3] * inverse;
          #pragma unroll
          for (BinaryInt j = 0; j<new_lin_size; ++j){
            const element_t el = matrix_list_old[pos4];
            matrix_list_new[pos1] = el;
            matrix_list_new[pos5] = el - matrix_list_old[pos6] * inverse_1;
            ++pos1;
            ++pos4;
            ++pos5;
            ++pos6;
          }
          ++pos4;
          pos3+=old_lin_size;
        }
      }
      for (BinaryInt set=0; set<two_level*2*new_mat_size; set+=new_mat_size){
        minors[mnr++] += matrix_list_new[set];
      }
      // print_matrix<2>(matrix_list_new, new_lin_size*(1<<(level+1)), new_lin_size);
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


}//namespace
