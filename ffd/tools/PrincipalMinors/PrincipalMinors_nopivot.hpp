namespace ffd::principal_minors{

  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  PrincipalMinors_nopivot(matrix_t&& matrix){
    

    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

    using namespace ffd::nilpotent_polynomial;
    using ffd::vector_range::Range;
    using std::abs;
    using std::cerr;
    using std::endl;
    using ffd::math_tools::print_matrix;

    const int linear_size = ffd::core_math::sqrt_int(size(matrix));
    assert (linear_size*linear_size == (int)size(matrix));

    const int cardinality = (1 << linear_size);

    std::vector<element_t> matrix_list((linear_size*(1+linear_size)*(1+2*linear_size)/6));
    

    std::vector<element_t> minors(cardinality);
    minors[0] = 1.0;
    minors[1] = matrix_list[0];

    // ffd::math_tools::print_matrix(matrix);

    // std::vector<BinaryInt> matrix_pos(order+1,0);
    // for (BinaryInt m_size=order-1; m_size>0; --m_size){
    //   matrix_pos[m_size] = matrix_pos[m_size+1] + (m_size+1)*(m_size+1);
    // }

    element_t inverse;
    element_t inverse_1;
    BinaryInt m_pos = 0;
    BinaryInt k_pos = linear_size - 1;
    int set = 1;
    std::vector<int> visited(linear_size+1,0U);
    visited[linear_size - 1] = 1;

    while(visited[0]!=cardinality-1)
    {
      if (k_pos==0) //bottom
      {
        k_pos++;
        set = visited[k_pos];
        BinaryInt k_pos_1 = k_pos+1;
        m_pos-=k_pos_1*k_pos_1;
      }
      else if (visited[k_pos-1]<set+(1<<(linear_size-k_pos-1))) //go left
      {
        BinaryInt k_pos_1 = k_pos+1;
        BinaryInt pos_ = m_pos+k_pos_1+1;
        m_pos += k_pos_1*k_pos_1;
        BinaryInt pos = m_pos;
        for (BinaryInt i = 0; i<k_pos; ++i){
          for (BinaryInt j = 0; j<k_pos; ++j){
            matrix_list[pos] = matrix_list[pos_];
            pos++;
            pos_++;
          }
          pos_++;
        }
        set += (1<<(linear_size-k_pos-1));
        k_pos--;
        visited[k_pos] = set;
        minors[set] = matrix_list[m_pos];
      }
      else if (visited[k_pos-1]<set+(1<<(linear_size-k_pos))) //go right
      {
        inverse = 1.0/minors[set];
        BinaryInt k_pos_1 = k_pos+1;
        BinaryInt pos_ = m_pos+k_pos_1+1;
        BinaryInt pos__ = m_pos+1;
        BinaryInt pos_i = m_pos+k_pos_1;
        m_pos += k_pos_1*k_pos_1;
        for (BinaryInt i = 0; i<k_pos; ++i){
          BinaryInt pos = m_pos+i;
          inverse_1 = matrix_list[pos_i] * inverse;
          for (BinaryInt j = 0; j<k_pos; ++j){
            matrix_list[pos] = matrix_list[pos_] - matrix_list[pos__] * inverse_1;
            pos+=k_pos;
            pos_++;
            pos__++;
          }
          pos_++;
          pos__-=k_pos;
          pos_i+=k_pos_1;
        }
        set += (1<<(linear_size-k_pos));
        k_pos--;
        visited[k_pos] = set;
        minors[set] = matrix_list[m_pos];
      }
      else //go up
      {
        visited[k_pos-1] = 0;
        k_pos++;
        set = visited[k_pos];
        BinaryInt k_pos_1 = k_pos+1;
        m_pos-= k_pos_1*k_pos_1;
      }
    }
    minors[cardinality-1] = matrix_list[m_pos];

    // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
    BinaryInt set_K = 1;
    for (int level = 0; level < linear_size; ++level){
      BinaryInt set_2K = 2*set_K;
      for (BinaryInt set = set_K+1; set < set_2K; ++set) {
        minors[set] *= minors[set-set_K];
      }
      set_K=set_2K;
    }

    // for (int set=0; set<(1<<linear_size);++set){
    //   std::cout << set << " " << minors[set] << std::endl;
    // }

    return minors;

  }

}//namespace
