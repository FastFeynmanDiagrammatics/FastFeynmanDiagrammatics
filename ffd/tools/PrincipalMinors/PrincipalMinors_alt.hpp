namespace ffd::principal_minors{


  // template<uint max_order = 0, typename matrix_t>
  template<typename matrix_t>
  std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
  PrincipalMinors(matrix_t&& matrix){


    using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

    using namespace ffd::nilpotent_polynomial;
    using ffd::vector_range::Range;
    using std::abs;
    using std::cerr;
    using std::endl;
    using ffd::math_tools::print_matrix;

    const int linear_size = ffd::core_math::sqrt_int(size(matrix));
    assert (linear_size*linear_size == (int)size(matrix));

    // uint max_ord;
    // if (max_order == 0){
    //   max_ord = linear_size;
    // }
    // else{
    //   max_ord = max_order;
    // }
    // cerr << "hre" << endl;

    const int cardinality = (1 << linear_size);

    std::vector<element_t> matrix_list((linear_size*(1+linear_size)*(1+2*linear_size)/6));

    std::vector<element_t> minors(cardinality);
    minors[0] = 1.0;

    // cerr << "begin" << endl;
    // for (int j : Range(cardinality)){
    //     for (int set : Range(minors[j].size())){
    //       cerr << j << " " << set << " " <<  minors[j][set] << endl;
    //     }
    // }
    // cerr << "end" << endl;

    element_t pivot = 0.0;
    for (BinaryInt element = 0; element < size(matrix); ++element){
      matrix_list[element] = matrix[element];
      pivot += abs(matrix[element]);
    }
    pivot /= (Real)size(matrix);
    // pivot = 0.0;
    // std::cerr << "pivot " << pivot << std::endl;

    // BinaryInt pop_set = 1;

    element_t inverse;
    element_t inverse_1;
    BinaryInt m_pos = 0;
    BinaryInt k_pos = linear_size - 1;
    int set = 1;
    std::vector<int> visited(linear_size+1,0U);
    // if constexpr(max_order!=0){
    //   visited[linear_size] = -1;
    // }
    visited[linear_size - 1] = 1;

    while(visited[0]!=cardinality-1)
    {
      // std::cerr <<"set" << set << endl;
      // if constexpr(max_order!=0){
      //   if (pop_set>max_order){
      //     minors[set] = matrix_list[m_pos] + pivot;
      //   }
      // }
      // else{
        minors[set] = matrix_list[m_pos] + pivot;
      // }

      // if constexpr(max_order!=0){ // bound
      //   if (pop_set>max_order){
      //     k_pos++;
      //     set = visited[k_pos];
      //     BinaryInt k_pos_1 = k_pos+1;
      //     m_pos-=k_pos_1*k_pos_1;
      //     pop_set--;
      //     continue;
      //   }
      // }

      if (k_pos==0) //bottom
      {
        k_pos++;
        set = visited[k_pos];
        BinaryInt k_pos_1 = k_pos+1;
        m_pos-=k_pos_1*k_pos_1;
        // if constexpr(max_order!=0){
        //   pop_set--;
        // }
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
        // if constexpr(max_order!=0){
        //   pop_set++;
        // }
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

    // if constexpr(max_order!=0){
    //   if (pop_set>max_order){
    //     minors[cardinality-1] = matrix_list[m_pos] + pivot;
    //   }
    // }
    // else{
      minors[cardinality-1] = matrix_list[m_pos] + pivot;
    // }

    // cerr << "hrse" << endl;

    //   cerr << "dfs" << endl;
    // for (BinaryInt set=0; set<cardinality; ++set){
    //   cerr << minors[set] << endl;
    // }

    // cerr << "begin" << endl;
    // for (int j : Range(cardinality)){
    //     for (int set : Range(minors[j].size())){
    //       cerr << j << " " << set << " " <<  minors[j][set] << endl;
    //     }
    // }
    // cerr << "end" << endl;

    // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
    BinaryInt set_K = 1;
    for (int level = 0; level < linear_size; ++level){
      BinaryInt set_2K = 2*set_K;
      for (BinaryInt set = set_K+1; set < set_2K; ++set) {
        minors[set] *= minors[set-set_K];
      }
      set_K=set_2K;
    }

    // PIVOT CORRECTIONS
    set_K = cardinality;
    for (int level = linear_size; level > 0; --level){
      int set_K_2 = (set_K>>1);
      for (int set = set_K-1; set >= set_K_2; --set){
        for (int set_ = set; set_ < cardinality; set_+=set_K){
          minors[set_] -= pivot * minors[set_-set_K_2];
        }
      }
      set_K = set_K_2;
    }

    // print_matrix<matrix_t>(minors,cardinality,1,0,"minors:");

    return minors;

  }

}//namespace
