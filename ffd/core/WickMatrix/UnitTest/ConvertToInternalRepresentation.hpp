

namespace ffd::wick_matrix::unit_test{

  void Convert_dimension(int dimension){
    using std::size;
    WickMatrix<int> M(Phi_(0), dimension);
    int component = 0;
    for(std::size_t j=0; j < size(M); ++j){
      for(std::size_t k=j; k < size(M); k++){
	M(j, k, 'a') = component;
	component++;
      }
    }
    for(std::size_t j=0; j < size(M.Components); ++j){
      assert(M.Components[j] == int(j));
    }

    WickMatrix<int> L(Rho_(0), dimension);
    component = 0;
    for(std::size_t j=0; j < size(L); ++j){
      for(std::size_t k=j+1; k < size(L); k++){
	L(j, k, 'a') = component;
	component++;
      }
    }
    for(std::size_t j=0; j < size(L.Components); ++j){
      assert(L.Components[j] == int(j) );
    }
  }

  
  void ConvertToInternalRepresentation(){
    Convert_dimension(2);
    Convert_dimension(3);
    Convert_dimension(4);
    Convert_dimension(10);
    Convert_dimension(11);
  }

}//namespace
