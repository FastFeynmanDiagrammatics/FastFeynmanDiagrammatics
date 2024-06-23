

namespace ffd::eigenvalues_vectors::unit_test{

  auto two_by_two_hermitian_matrix_test_arg(std::array<std::array<Complex, 2>, 2> H_){
    bool IsOk = true;

    
    
    using std::real, std::imag, std::sqrt, std::abs, std::atan2, std::asin;

    
    
    Real phi = std::arg(H_[0][1]);
    auto P_phi = PhaseMatrix_2x2(phi);
    auto P_phi_dagger = PhaseMatrix_2x2(-phi);
    auto H_rotated = MultiplyMatrices_array(H_, P_phi);
    H_rotated = MultiplyMatrices_array(P_phi_dagger, H_rotated);


    Real epsilon = numeric_limits<Real>::epsilon();
    for(int j: {0, 1}){
      for(int k: {0, 1}){
	IsOk = IsOk && abs( imag( H_rotated[j][k] ) ) < 10*epsilon;
	IsOk = IsOk && abs( H_rotated[j][k] - H_rotated[k][j]) < 10*epsilon;
      }
    }
    
    return IsOk;
  }


}//namespace
