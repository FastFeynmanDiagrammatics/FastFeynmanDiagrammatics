#ifndef FFD_PFAFFIAN_HEADER_NOT_BEEN
#define FFD_PFAFFIAN_HEADER_NOT_BEEN


#include"../Matrix.hpp"


namespace FFD::Pfaffian{

  template<typename Field>
  using Matrix = FFD::Matrix::Matrix<Field>;
  
  template<typename Field>
  Field Pfaffian(Matrix<Field> M){
    if(M.Size()==0){
      return 1;
    }
    if(M.Size()%2 == 1){
      return 0;
    }
    return 1; 
  }

  
}//namespace
#endif
