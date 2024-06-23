

namespace ffd::eigenvalues_vectors::unit_test{

  template<typename Field>
  bool eigenvalues_vectors_matrix_assertable(std::vector<std::vector<Field>> matrix,
					     Real precision = -1.,
					     bool DontPrint = true,
					     Real precision_check = -1.){
			
    bool IsOk = true;

    
    Real Precision = precision > 0 ?
      precision :
      20*std::numeric_limits<Real>::epsilon() ;


    Real PrecisionCheck = precision_check > 0 ?
      precision_check :
      Precision;
    
    
    auto [Eigenvalues, Eigenvectors] =
      EigenvaluesVectors( matrix, Precision );
    int const size = std::size(Eigenvalues);
    Real Norm2_matrix = MatrixNorm<2>(matrix);


    for(int m = 0; m < size; ++m){
      auto Eigenvector = Eigenvectors[m];

      
      Real Norm2_Eigenvector = 0.;
      for(int k = 0; k < size; ++k){
	Norm2_Eigenvector += pow(std::abs(Eigenvector[k]), 2);
      }
      if(!DontPrint){
	std::cerr<<std::setprecision(std::numeric_limits<Real>::digits10 );
	std::cerr<<"norm="<<Norm2_Eigenvector<<std::endl;
      }
      IsOk = IsOk &&
	std::abs(Norm2_Eigenvector-1.) <
	PrecisionCheck;


      auto lambda = Eigenvalues[m];
      if(!DontPrint){
	std::cerr<<"lambda="<<lambda<<std::endl;
      }
      
      std::vector<Field> lambda_Eigenvector(size, 0.);
      for(int j = 0; j < size; ++j){
	for(int k = 0; k < size; ++k){
	  lambda_Eigenvector[j] += matrix[j][k]*Eigenvector[k];
	}
	if(!DontPrint){
	  std::cerr<<lambda_Eigenvector[j]<<" "<<lambda*Eigenvector[j]-lambda_Eigenvector[j]<<std::endl;
	}
	IsOk = IsOk &&
	  std::abs(lambda_Eigenvector[j] - lambda*Eigenvector[j]) <
	  Norm2_matrix*PrecisionCheck;
      }
    }
    return IsOk;
  }


}//namespace
