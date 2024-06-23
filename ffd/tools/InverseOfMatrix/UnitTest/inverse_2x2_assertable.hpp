

namespace ffd::inverse_matrix::unit_test{

  template<typename Field>
  auto
  inverse_2x2_assertable(std::vector<Field> matrix,
			 Real precision = -1.){
    bool IsOk = true;

    

    Real Precision = precision > 0 ?
      precision :
      10*std::numeric_limits<Real>::epsilon();

    

    assert( std::size(matrix) == 4 );

    

    decltype(matrix) inverse(4, 0.);

    

    Field det = matrix[0]*matrix[3]-matrix[1]*matrix[2];

    

    for(int j: {0, 1}){
      inverse[j+2*j] =      matrix[1-j+2*(1-j)]/det;
      inverse[j+2*(1-j)] = -matrix[j+2*(1-j)]  /det;
    }

    


    auto inverse2 = InverseOfMatrix<Field>(matrix);


    
    for(int j: {0, 1, 2, 3}){
      IsOk = IsOk &&
	std::abs(inverse[j] - inverse2[j]) <
	Precision;
    }
    return IsOk;
  }
  


}//namespace
