

namespace ffd::principal_minors_old{

template<typename Field>
std::vector<Field>
CalculateMinors(std::vector<Field> matrix){

  using namespace ffd::nilpotent_polynomial;
  using ffd::vector_range::Range;
  using std::abs;
  using std::cerr;
  using std::endl;

  const int linear_size = ffd::core_math::sqrt_int(size(matrix));
  assert (linear_size*linear_size == (int)size(matrix));
\
  const int cardinality = (1 << linear_size);

  std::vector<Field> matrix_list((linear_size*(1+linear_size)*(1+2*linear_size)/6));

  std::vector<Field> minors(cardinality);
  minors[0] = 1.0;

  // cerr << "begin" << endl;
  // for (int j : Range(cardinality)){
  //     for (int set : Range(minors[j].size())){
  //       cerr << j << " " << set << " " <<  minors[j][set] << endl;
  //     }
  // }
  // cerr << "end" << endl;

  Field pivot = 0.0;
  for (unsigned int element = 0; element < size(matrix); ++element){
    matrix_list[element] = matrix[element];

    if constexpr(std::is_same_v<Real,Field> || std::is_same_v<Complex,Field>)
    {
      pivot += abs(matrix[element]);
    }
    else if constexpr(std::is_same_v<NilpotentPolynomial<Real>,Field > || std::is_same_v<NilpotentPolynomial<Complex>,Field>){
      pivot += poly_abs(matrix[element]);
      // pivot += abs(matrix[element]);
    }
    else{
      assert(false);
    }

  }
  pivot /= (Real)size(matrix);
  // std::cerr << "pivot " << pivot << std::endl;

  Field inverse;
  Field inverse_1;
  unsigned int m_pos = 0;
  unsigned int k_pos = linear_size - 1;
  unsigned int set = 1;
  std::vector<unsigned int> visited(linear_size,0U);
  visited[linear_size - 1] = 1;

    while((int)visited[0]!=cardinality-1)
    {
      minors[set] = matrix_list[m_pos] + pivot;
      if (k_pos==0) //bottom
      {
        k_pos++;
        set = visited[k_pos];
        unsigned int k_pos_1 = k_pos+1;
        m_pos-=k_pos_1*k_pos_1;
      }
      else if (visited[k_pos-1]<set+(1<<(linear_size-k_pos-1))) //go left
      {
        unsigned int k_pos_1 = k_pos+1;
        unsigned int pos_ = m_pos+k_pos_1+1;
        m_pos += k_pos_1*k_pos_1;
        unsigned int pos = m_pos;
        for (unsigned int i = 0; i<k_pos; ++i)
        {
          for (unsigned int j = 0; j<k_pos; ++j)
          {
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
        inverse = 1.0/minors[set]; //?
        unsigned int k_pos_1 = k_pos+1;
        unsigned int pos_ = m_pos+k_pos_1+1;
        unsigned int pos__ = m_pos+1;
        unsigned int pos_i = m_pos+k_pos_1;
        m_pos += k_pos_1*k_pos_1;
        for (unsigned int i = 0; i<k_pos; ++i)
        {
          unsigned int pos = m_pos+i;
          inverse_1 = matrix_list[pos_i] * inverse; //?
          for (unsigned int j = 0; j<k_pos; ++j)
          {
            matrix_list[pos] = matrix_list[pos_] - matrix_list[pos__] * inverse_1; //?
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
      }
      else //go up
      {
        visited[k_pos-1] = 0;
        k_pos++;
        set = visited[k_pos];
        unsigned int k_pos_1 = k_pos+1;
        m_pos-= k_pos_1*k_pos_1;
      }
    }

    minors[cardinality-1] = matrix_list[m_pos] + pivot;

    // cerr << "begin" << endl;
    // for (int j : Range(cardinality)){
    //     for (int set : Range(minors[j].size())){
    //       cerr << j << " " << set << " " <<  minors[j][set] << endl;
    //     }
    // }
    // cerr << "end" << endl;

    // MULTIPLY BY ELEMENT OF PREVIOUS LEVEL
    unsigned int set_K = 1;
    for (int level = 0; level < linear_size; ++level)
    {
      unsigned int set_2K = 2*set_K;
      for (unsigned int set = set_K+1; set < set_2K; ++set) // changed this recently from set = set_K
      {
        unsigned int alt_set = set-set_K;
        minors[set] *= minors[alt_set];
      }
      set_K=set_2K;
    }

    // PIVOT CORRECTIONS
    set_K = cardinality;
    for (int level = linear_size; level > 0; --level)
    {
      int set_K_2 = (set_K>>1);
      for (int set = set_K-1; set >= set_K_2; --set)
      {
        for (int set_ = set; set_ < cardinality; set_+=set_K)
        {
          int alt_set = set_-set_K_2;
          minors[set_] -= pivot * minors[alt_set];
        }
      }
      set_K = set_K_2;
    }

  return minors;

}

}//namespace
