#ifndef FFD_MATRIX_HEADER_NOT_BEEN_HERE
#define FFD_MATRIX_HEADER_NOT_BEEN_HERE


namespace ffd::Matrix{

  template<int degree, typename Field, int nilpotency>
  using NilpotentPolynomial = ffd::NilpotentPolynomial::NilpotentPolynomial<degree, Field, nilpotency>;

  
  template <typename Field>
  Field Pivot(const Field& x){
    return x;
  }


  template <int degree, typename Field, int nilpotency>
  Field Pivot(const NilpotentPolynomial<degree, Field, nilpotency>& P1){
    return P1.coef[0];
  }
  
  
  enum class MatrixSymmetry{NoSymmetry, Antisymmetric, Symmetric};


  template<typename Field>
  struct Matrix{
    std::vector<Field> coef;
    int Size;
    MatrixSymmetry Symmetry;
    
    Matrix(int size=0, MatrixSymmetry symmetry = MatrixSymmetry::NoSymmetry): Size(size), Symmetry(symmetry),
      coef ((symmetry == MatrixSymmetry::NoSymmetry)*size*size +
      (symmetry == MatrixSymmetry::Antisymmetric)*(size*(size-1))/2 +
	    (symmetry == MatrixSymmetry::Symmetric)*(size*(size+1))/2, 0){}
    
  
    Field& operator()(int j, int k){
      if(j >= k && Symmetry == MatrixSymmetry::Antisymmetric){
	throw std::invalid_argument("Matrix::Matrix::operator(): #1 >= #2 for antisymmetric matrices is not supported");
      }
      if(Symmetry == MatrixSymmetry::NoSymmetry){
	return coef[k+j*Size];
      }else if(Symmetry == MatrixSymmetry::Antisymmetric){
	int previous_rows = 0;
	if(j>0){
	  previous_rows += Size*j - (j*(j+1))/2;
	}
	int k_delta = k-j-1;
	return coef[previous_rows + k_delta];
      }else{
	if(j>k){
	  std::swap(j,k);
	}
	int previous_rows = 0;
	if(j>0){
	  previous_rows += Size*j - (j*(j-1))/2;
	}
	int k_delta = k-j;
	return coef[previous_rows + k_delta];
      }
    }
  };

  
}//namespace FFD::Matrix




#endif
