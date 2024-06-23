

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<int d, int j, typename T>
  void
  SetSpaceCoordinate(T& X, std::array<int, d> r){
    component<j+1>(X) = r[j];


    if constexpr( j+1 < d){
	SetSpaceCoordinate<d, j+1>(X, r);
      }
  }

  
  template<int d, typename T>
  void
  SetSpaceCoordinates(T& X, std::array<int, d> r){
    SetSpaceCoordinate<d, 0>(X, r);
  }
  
  

}//namespace
