

namespace ffd::eigenvalues_vectors::unit_test{

  auto RotationMatrix_2x2(Real theta){
    std::array<std::array<Real, 2>, 2> rotation_matrix;
    
    for(int j: {0, 1}){
      rotation_matrix[j][j] = std::cos(theta);
      rotation_matrix[j][1-j] = (2*j-1)*std::sin(theta);
    }
    return rotation_matrix;
  }


}//namespace
