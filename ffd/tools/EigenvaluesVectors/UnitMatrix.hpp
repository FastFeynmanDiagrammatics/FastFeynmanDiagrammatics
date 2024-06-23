

namespace ffd::eigenvalues_vectors{

  template<typename Field>
  auto UnitMatrix(int size){
    std::vector<std::vector<Field>>
      unit_matrix(size, std::vector<Field>(size, 0.));

    
    for(int j = 0; j < size; ++j){
      unit_matrix[j][j] = 1.;
    }
    return unit_matrix;
  }
  

}//namespace
