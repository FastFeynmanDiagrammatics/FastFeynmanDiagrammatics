

namespace ffd::eigenvalues_vectors::unit_test{

  template<typename Field>
  bool inverse_assertable_2x2(std::vector<std::vector<Field>> M,
			      Real precision = -1){
    bool IsOk = true;



    Real Precision = precision > 0 ?
      precision :
      10*std::numeric_limits<Real>::epsilon();

    

    auto Inverse = InverseMatrix(M);

    
    
    std::vector<std::vector<Field>> inverse(2, std::vector<Field>(2, 0.));

    
    
    Field det = M[0][0]*M[1][1] - M[0][1]*M[1][0];


    
    for(int j: {0, 1}){
      inverse[j][j]   = M[1-j][1-j]/det;
      inverse[j][1-j] = -M[j][1-j]/det;
    }

    for(int j: {0, 1}){
      for(int k: {0, 1}){
	IsOk = IsOk &&
	  std::abs( inverse[j][k] - Inverse[j][k] )
	  < Precision;
      }
    }
    return IsOk;
  }
  


}//namespace
