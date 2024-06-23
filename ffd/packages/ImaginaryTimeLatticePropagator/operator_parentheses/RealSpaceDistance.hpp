

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<int d, typename T>
  
  auto

  RealSpaceDistance(ffd::feynman_edge::
		    QuantumFieldPositions const& Psi_X){
    std::array<std::array<Real, d>, 2> real_spaces;

    
    using ffd::vector_range::Range;
    

    for( int j: {0, 1} ){
      real_spaces[j] =
	ffd::RealSpace(
		       std::any_cast<T>(
					std::get<std::any>( Psi_X[j] )
					)
		       );
    }
    
    
    std::array<Real, d> real_space_distance;
    real_space_distance.fill(0.);
    for(  auto [j, u]: ffd::user_space::
	    CartesianProductRange( std::array<int, 2>{d, 2} )  ){
      real_space_distance[j] += (1-2*u)*real_spaces[u][j];
    }


    return real_space_distance;
  }

}//namespace
