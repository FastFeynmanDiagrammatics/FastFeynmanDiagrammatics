

namespace ffd::eigenvalues_vectors{

  template<typename Field>
  auto
  InverseMatrix(std::vector<std::vector<Field>> const& M_){
    int const size = std::size(M_);
    std::vector<std::vector<Field>> Inverse(size, std::vector<Field>(size, 0.) );
    
    

    
    std::vector<std::vector<Complex>> U(size, std::vector<Complex>(size, 0.) );
    

    
    auto [eigenvalues, eigenvectors] =
      EigenvaluesVectors(M_);



    for(int j=0; j < size; ++j){
      for(int k=0; k < size; ++k){
	U[j][k] = std::conj(eigenvectors[j][k]);
	// U[j][k] = eigenvectors[j][k];
      }
    }

    

    std::vector<Complex> complex_diagonal(size, 0.);
    for(int j=0; j< size; ++j){
      complex_diagonal[j] = 1./eigenvalues[j];
    }

    

    auto diagonal_matrix = DiagonalMatrix(complex_diagonal);


    
    auto Complex_Inverse = MultiplyMatrices(diagonal_matrix,
					    U);

        

    Complex_Inverse = MultiplyMatrices(Dagger(U),
				       Complex_Inverse);

    


    if constexpr(std::is_same<Field, Real>::value ){
	for(int j=0; j < size; ++j){
	  for(int k=0; k < size; ++k){
	    Inverse[j][k] = std::real(Complex_Inverse[j][k]);
	  }
	}
      }else{
      Inverse = Complex_Inverse;
    }
    return Inverse;
  }


}//namespace
