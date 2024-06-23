

namespace ffd::wick_matrix::unit_test{
  
  void OperatorParentheses(){
    using std::size;
    WickMatrix<int> M(Psi_(1), 10);
    for(std::size_t j=0; j < size(M); ++j){
      for(std::size_t k=0; k < size(M); ++k){
	M(j, k, 'a') = k+j*size(M);
	assert( M(j, k) == int(k+j*size(M)) );
      }
    }

    WickMatrix<int> M2(Eta_(1), 10);
    for(std::size_t j=0; j < size(M2); ++j){
      for(std::size_t k=0; k < size(M2); ++k){
	M2(j, k, 'a') = k+j*size(M2);
	assert( M2(j, k) == int(k+j*size(M2)) );
      }
    }

    WickMatrix<int> M3(Phi_(1), 10);
    for(std::size_t j=0; j < size(M3); ++j){
      for(std::size_t k=0; k < size(M3); ++k){
	M3(j, k, 'a') = k+j*size(M3);
	assert( M3(j, k) == int(k+j*size(M3)) );
	assert( M3(j, k) == M3(k, j) );
      }
    }

    WickMatrix<int> M4(Rho_(1), 10);
    for(std::size_t j=0; j < size(M4); ++j){
      for(std::size_t k=j+1; k < size(M4); ++k){
	M4(j, k, 'a') = k+j*size(M4);
	assert( M4(j, k) == int(k+j*size(M4)) );
	assert( M4(j, k) == - M4(k, j) );
      }
    }
  }

}//namespace
