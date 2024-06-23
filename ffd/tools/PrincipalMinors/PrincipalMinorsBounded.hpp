namespace ffd::principal_minors{

  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  PrincipalMinorsBounded(matrix_t&& matrix, const BinaryInt max_order){


    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

    using namespace ffd::nilpotent_polynomial;
    using ffd::vector_range::Range;
    using std::abs;
    using std::cerr;
    using std::cout;
    using std::endl;
    using ffd::math_tools::print_matrix;

    const int order = ffd::core_math::sqrt_int(size(matrix));
    assert (order*order == (int)size(matrix));

    const int cardinality = (1 << order);

    std::vector<element_t> matrix_list((order*(1+order)*(1+2*order)/6));

    element_t pivot = 0.0;
    for (BinaryInt element = 0; element < size(matrix); ++element){
      matrix_list[element] = matrix[element];
      pivot += abs(matrix[element]);
    }
    pivot /= (Real)size(matrix);

    std::vector<element_t> minors(cardinality);
    minors[0] = 1.0;
    minors[1] = matrix_list[0]+pivot;

    element_t inverse;
    element_t inverse_1;
    BinaryInt m_pos = 0;
    BinaryInt k_pos = order - 1;
    int set = 1;
    std::vector<int> visited(order+1,0U);
    visited[order - 1] = 1;
    BinaryInt pop_set = 1;

    // BinaryInt max_set = (1<<order) - (1<<(order-max_order));

    while(pop_set>0)
    {
      if (k_pos==0) //bottom
      {
        k_pos++;
        set = visited[k_pos];
        BinaryInt k_pos_1 = k_pos+1;
        m_pos-=k_pos_1*k_pos_1;
      }
      else if (pop_set>max_order) //bound
      {
        // cout << "bound " << set << endl;
        minors[set] = 0;
        visited[k_pos-1] = 0;
        k_pos++;
        set = visited[k_pos];
        BinaryInt k_pos_1 = k_pos+1;
        m_pos-= k_pos_1*k_pos_1;
      }
      else if (visited[k_pos-1]<set+(1<<(order-k_pos-1))) //go left
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
        set += (1<<(order-k_pos-1));
        k_pos--;
        visited[k_pos] = set;
        minors[set] = matrix_list[m_pos] + pivot;
      }
      else if (visited[k_pos-1]<set+(1<<(order-k_pos))) //go right
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
        set += (1<<(order-k_pos));
        k_pos--;
        visited[k_pos] = set;
        if (__builtin_popcount(set)<=max_order)
        minors[set] = matrix_list[m_pos] + pivot;
        // pop_set++;
      }
      else //go up
      {
        visited[k_pos-1] = 0;
        k_pos++;
        set = visited[k_pos];
        BinaryInt k_pos_1 = k_pos+1;
        m_pos-= k_pos_1*k_pos_1;
      }
      pop_set = __builtin_popcount(set);
    }

    // for (BinaryInt set =0; set< cardinality; ++set) {
    //   cout << set << " " << minors[set] << endl;
    // }

    // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
    BinaryInt set_K = 1;
    for (int level = 0; level < order; ++level){
      BinaryInt set_2K = 2*set_K;
      for (BinaryInt set = set_K+1; set < set_2K; ++set) {
        minors[set] *= minors[set-set_K];
      }
      set_K=set_2K;
    }

    // PIVOT CORRECTIONS
    set_K = cardinality;
    for (int level = order; level > 0; --level){
      int set_K_2 = (set_K>>1);
      for (int set = set_K-1; set >= set_K_2; --set){
        for (int set_ = set; set_ < cardinality; set_+=set_K){
          if (__builtin_popcount(set_)<=max_order)
          minors[set_] -= pivot * minors[set_-set_K_2];
        }
      }
      set_K = set_K_2;
    }

    return minors;

  }

}//namespace
