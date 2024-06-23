

namespace ffd::eigenvalues_vectors{


  template<typename Field>
  auto
  RotationMatrix(Real theta,
		 std::array<int, 2> components_to_be_rotated,
		 Real size){
    std::vector<std::vector<Field>> rotation_matrix =
      UnitMatrix<Field>(size);


    
    auto const& k = components_to_be_rotated;
    assert(k[0] != k[1]);


    
    
    for(int j: {0, 1}){
      rotation_matrix[ k[j] ][ k[j]   ]  = std::cos(theta);
      rotation_matrix[ k[j] ][ k[1-j] ]  = (2*j-1)*std::sin(theta);
    }
    return rotation_matrix;
  }
  

}//namespace
