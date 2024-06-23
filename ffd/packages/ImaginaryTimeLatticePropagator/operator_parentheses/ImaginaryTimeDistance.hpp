

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<typename T>
  
  auto
  
  ImaginaryTimeDistance(ffd::feynman_edge::
			QuantumFieldPositions Psi_X){
    
    
    std::array<Real, 2> imag_times;
    for( int j: {0, 1} ){
     imag_times[j] =
       ffd::get<0>(
		   std::any_cast<T>(
				    std::get<std::any>(
						       Psi_X[j]
						       )
				    )
		   )();
    }
    
    
    Real diff_times = 0;
    for( int j: {0, 1} ){
      diff_times += (1-2*j)*imag_times[j];
    }


    return diff_times;

    
  }

}//namespace
