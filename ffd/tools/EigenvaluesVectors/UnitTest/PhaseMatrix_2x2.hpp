

namespace ffd::eigenvalues_vectors::unit_test{

  auto PhaseMatrix_2x2(Real phi_){
    std::array<std::array<Complex, 2>, 2> phase_matrix;

    
    using ffd::user_space::expI;


    for(int j: {0, 1}){
      phase_matrix[j][j] = expI((.5-j)*phi_);
      phase_matrix[j][1-j] = 0.;
    }
    return phase_matrix;
  }
		   

}//namespace
