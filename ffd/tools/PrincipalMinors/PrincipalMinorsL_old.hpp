namespace ffd::principal_minors{

template<typename matrix_t>
std::vector< typename std::decay<decltype(std::declval<matrix_t>()[0])>::type >
PMOLD(matrix_t&& matrix, BinaryInt step, Real epsilon = 10000*std::numeric_limits<Real>::min()){

  using element_t = typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
  using namespace ffd::nilpotent_polynomial;
  using ffd::vector_range::Range;
  using std::abs;
  // using std::cerr;
  // using std::endl;
  using ffd::math_tools::print_matrix;
  using ffd::inverse_matrix_gauss::Inverse_and_Determinant;

  const int order = ffd::core_math::sqrt_int(size(matrix)/(step*step));
  assert (order*order*step*step == (int)size(matrix));
  const int card = (1 << order);

  std::vector<element_t> matrix_list(step*step*(order*(1+order)*(1+2*order)/6));
  std::vector<element_t> minors(card);
  minors[0] = 1.0;

  for (BinaryInt element = 0; element < size(matrix); ++element){
    matrix_list[element] = matrix[element];
  }

  std::vector<BinaryInt> matrix_pos(order+1,0);
  for (std::size_t m_size=order-1; m_size>0; --m_size){
	  matrix_pos[m_size] = matrix_pos[m_size+1] + (m_size+1)*(m_size+1)*step*step;
  }

  BinaryInt pos = matrix_pos[order];
  BinaryInt pos_1 = matrix_pos[order-1];

  BinaryInt level = order - 1;
  BinaryInt set = 1;
  std::vector<BinaryInt> visited(order,0U);
  visited[order - 1] = 1;

  std::vector<element_t> block(step*step,0.);

  // could replace at some point by: for(BinaryInt step=0; step<(1<<(order+1))-order-3; ++step)
  while(visited[0]!=card-1){
	  		// std::cerr << set << " old set" << std::endl;
    pos = matrix_pos[level+1];
    for (BinaryInt col=0; col<step; ++col){
      for (BinaryInt row=0; row<step; ++row){
        block[col*step+row]= matrix_list[pos+col*(level+1)*step+row];
      }
    }

    auto [inverse,determinant] = Inverse_and_Determinant(block);
    if (abs(determinant) < epsilon) return DeterminantPrincipalMinors(matrix, step);

    minors[set] = determinant;
// std::cerr << set << " old minor " << minors[set] << std::endl;
    // BOTTOM
    if (level==0){
      level++;
      set = visited[level];
    }
    // GO LEFT
    else if (visited[level-1]<set+(1<<(order-level-1))){

  		pos = matrix_pos[level+1] + (level+1)*step*step + step;
  		pos_1 = matrix_pos[level];
      for (BinaryInt i = 0; i<level*step; ++i){
        for (BinaryInt j = 0; j<level*step; ++j){
          matrix_list[pos_1] = matrix_list[pos];
          pos++;
          pos_1++;
        }
        pos+=step;
      }

      set += (1<<(order-level-1));
      level--;
      visited[level] = set;
    }
    // GO RIGHT
    else if (visited[level-1]<set+(1<<(order-level))){

      std::vector<element_t> temp_1(step*step*level,0.);
      std::vector<element_t> temp_2(step*step*level*level,0.);

      for (BinaryInt i=0; i<step; ++i){
    		pos = matrix_pos[level+1]+step;
        for (BinaryInt k=0; k<step; ++k){
          for (BinaryInt j=0; j<step*level; ++j){
            temp_1[i*step*level+j] += inverse[i*step+k] * matrix_list[pos];
            pos++;
          }
          pos+=step;
        }
      }

  		pos = matrix_pos[level+1]+step*step*(level+1);
      for (BinaryInt i=0; i<step*level; ++i){
        for (BinaryInt k=0; k<step; ++k){
          for (BinaryInt j=0; j<step*level; ++j){
            temp_2[i*step*level+j] += temp_1[k*level*step+j] * matrix_list[pos];
          }
          pos++;
        }
        pos+=step*level;
      }

      pos = matrix_pos[level+1]+step*step*(level+1)+step;
      pos_1 = matrix_pos[level];
      for (BinaryInt i=0; i<level*step; ++i){
        for (BinaryInt j=0; j<level*step; ++j){
          matrix_list[pos_1 + i*level*step+j] = matrix_list[pos + i*level*step+j] - temp_2[i*step*level+j];
        }
        pos+=step;
      }

      set += (1<<(order-level));
      level--;
      visited[level] = set;
    }
    // GO UP
    else{
      visited[level-1] = 0;
      level++;
      set = visited[level];
    }
  }
  pos = matrix_pos[level+1];
    // std::cerr << "old " << pos << " " << level << " " << std::endl;
  for (BinaryInt col=0; col<step; ++col){
    for (BinaryInt row=0; row<step; ++row){
      block[col*step+row]= matrix_list[pos+col*(level+1)*step+row];
    }
  }
  auto [inverse,determinant] = Inverse_and_Determinant(block);
  if (abs(determinant) < epsilon) return DeterminantPrincipalMinors(matrix, step);
  minors[card-1] = determinant;
// std::cerr << card-1 << " old minor " << minors[card-1] << std::endl;

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
//  std::cerr <<"ret "<<std::endl;
  return minors;
}

}//namespace
