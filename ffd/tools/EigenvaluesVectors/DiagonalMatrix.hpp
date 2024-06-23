

namespace ffd::eigenvalues_vectors{

  template<typename Field>
  auto
  DiagonalMatrix(std::vector<Field> const& diagonal){
    int const size = std::size(diagonal);
    auto diagonal_matrix = UnitMatrix<Field>(size);


    for(int j=0; j < size; ++j){
      diagonal_matrix[j][j] = diagonal[j];
    }
    return diagonal_matrix;
  }

}//namespace
