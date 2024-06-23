namespace ffd::principal_minors{

  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  PrincipalMinors2n(matrix_t&& matrix, Real epsilon = 10000*std::numeric_limits<Real>::min()){

    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
    using namespace ffd::nilpotent_polynomial;
    using ffd::vector_range::Range;
    using std::abs;
    // using std::cerr;
    // using std::endl;
    // using ffd::math_tools::print_matrix;
    // using ffd::inverse_matrix_gauss::Inverse_and_Determinant;
    // using ffd::determinant::Determinant;

    const int order = ffd::core_math::sqrt_int(size(matrix)/(4));
    assert (order*order*4 == (int)size(matrix));
    const int card = (1 << order);

    std::vector<element_t> matrix_list(4*(order*(1+order)*(1+2*order)/6));
    std::vector<element_t> minors(card);
    minors[0] = 1.0;

    for (BinaryInt element = 0; element < size(matrix); ++element){
      matrix_list[element] = matrix[element];
    }

    std::vector<BinaryInt> matrix_pos(order+1,0);
    for (std::size_t m_size=order-1; m_size>0; --m_size){
      matrix_pos[m_size] = matrix_pos[m_size+1] + (m_size+1)*(m_size+1)*4;
    }

    std::vector<element_t> inverse(4,0.);
    BinaryInt pos = matrix_pos[order];
    BinaryInt pos_1 = matrix_pos[order-1];
    BinaryInt level = order - 1;
    BinaryInt set = 1;
    std::vector<BinaryInt> visited(order,0U);
    visited[order - 1] = 1;

    int l_1_2 = (level+1)*2;
    element_t determinant = matrix_list[0]*matrix_list[1+l_1_2] - matrix_list[1]*matrix_list[l_1_2];
    // if (abs(determinant) < epsilon) return DeterminantPrincipalMinors(matrix, 2);
    minors[set] = determinant;
    // std::cerr << set << " minor " << minors[set] << std::endl;

    // could replace at some point by: for(BinaryInt 2=0; 2<(1<<(order+1))-order-3; ++2)
    while(visited[0]!=card-1){
      // std::cerr << set << " set" << std::endl;

      // BOTTOM
      if (level==0){

        level++;
        set = visited[level];
        pos = matrix_pos[level+1];

      }

      // GO LEFT
      else if (visited[level-1]<set+(1<<(order-level-1))){

        BinaryInt l_1 = level+1;
        BinaryInt l_2 = 2*level;
        pos = matrix_pos[l_1] + 2*l_2 + 6;
        pos_1 = matrix_pos[level];
        for (BinaryInt i = 0; i<l_2; ++i){
          for (BinaryInt j = 0; j<l_2; ++j){
            matrix_list[pos_1] = matrix_list[pos];
            pos++;
            pos_1++;
          }
          pos+=2;
        }

        set += (1<<(order-level-1));
        level--;
        visited[level] = set;
        pos = matrix_pos[level+1];

        l_1_2 = (level+1)*2;
        determinant = matrix_list[pos+0]*matrix_list[pos+1+l_1_2]
        - matrix_list[pos+1]*matrix_list[pos+0+l_1_2];
        minors[set] = determinant;
        // std::cerr << set << " minor " << minors[set] << std::endl;

      }

      // GO RIGHT
      else if (visited[level-1]<set+(1<<(order-level))){

        BinaryInt s_l = 2*level;
        BinaryInt s_l_1 = s_l+2;

        inverse[0] = matrix_list[pos+1+s_l_1];  // d
        inverse[1] = -matrix_list[pos+1];       //-c
        inverse[2] = -matrix_list[pos+s_l_1];   //-b
        inverse[3] = matrix_list[pos+0];        // a
        determinant = inverse[0]*inverse[3] - inverse[1]*inverse[2];
        // if (abs(determinant) < epsilon) return DeterminantPrincipalMinors(matrix, 2);
        // if (abs(determinant) < epsilon) std::cout << "det issue" << std::endl;
        inverse[0] /= determinant;
        inverse[1] /= determinant;
        inverse[2] /= determinant;
        inverse[3] /= determinant;

        std::vector<element_t> temp_1(4*level,0.);
        std::vector<element_t> temp_2(4*level*level,0.);

        pos = matrix_pos[level+1]+2;
        for (BinaryInt j=0; j<s_l; ++j){
          element_t m1 = matrix_list[pos+j];
          element_t m2 = matrix_list[pos+s_l_1+j];
          temp_1[j]     += inverse[0] * m1 + inverse[1] * m2;
          temp_1[s_l+j] += inverse[2] * m1 + inverse[3] * m2;
        }

        pos = matrix_pos[level+1]+4*(level+1);
        for (BinaryInt i=0; i<s_l; ++i){
          BinaryInt i_s_l = i*s_l;
          element_t m_pos_0 = matrix_list[pos];
          element_t m_pos_1 = matrix_list[pos+1];
          for (BinaryInt j=0; j<s_l; ++j){
            temp_2[i_s_l+j] += temp_1[j] * m_pos_0 + temp_1[s_l+j] * m_pos_1;
          }
          pos+=s_l_1;
        }

        pos = matrix_pos[level+1]+4*(level+1)+2;
        pos_1 = matrix_pos[level];
        for (BinaryInt i=0; i<s_l; ++i){
          BinaryInt i_s_l = i*s_l;
          BinaryInt p_i_s_l = pos + i_s_l;
          BinaryInt p1_i_s_l = pos_1 + i_s_l;
          for (BinaryInt j=0; j<s_l; ++j){
            matrix_list[p1_i_s_l + j] = matrix_list[p_i_s_l + j] - temp_2[i_s_l + j];
          }
          pos+=2;
        }

        set += (1<<(order-level));
        level--;
        visited[level] = set;
        pos = matrix_pos[level+1];

        l_1_2 = (level+1)*2;
        determinant = matrix_list[pos+0]*matrix_list[pos+1+l_1_2]
        - matrix_list[pos+1]*matrix_list[pos+l_1_2];
        minors[set] = determinant;
        // std::cerr << set << " minor " << minors[set] << std::endl;

      }

      // GO UP
      else{

        visited[level-1] = 0;
        level++;
        set = visited[level];
        pos = matrix_pos[level+1];

      }
    }

    pos = matrix_pos[level+1];

    l_1_2 = (level+1)*2;
    determinant = matrix_list[pos+0]*matrix_list[pos+1+l_1_2]
    - matrix_list[pos+1]*matrix_list[pos+l_1_2];
    // if (abs(determinant) < epsilon) return DeterminantPrincipalMinors(matrix, 2);
    minors[card-1] = determinant;

    // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
    BinaryInt set_K = 1;
    for (int level = 0; level < order; ++level)
    {
      BinaryInt set_2K = 2*set_K;
      for (BinaryInt set = set_K+1; set < set_2K; ++set)
      {
        BinaryInt alt_set = set-set_K;
        minors[set] *= minors[alt_set];
      }
      set_K=set_2K;
    }

    return minors;
  }



}//namespace
