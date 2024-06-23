

namespace ffd::eigenvalues_vectors::unit_test{

  auto two_by_two_matrix_test_arg(std::array<std::array<Real, 2>, 2> S_){
    bool IsOk = true;

    
    using ffd::core_math::Pi, std::abs;
    using namespace std;
    
    
    Real sum = S_[0][0] + S_[1][1];
    Real diff = std::hypot(S_[0][0]-S_[1][1], 2*S_[0][1]);
    Real determinant = S_[0][0]*S_[1][1] - S_[0][1]*S_[1][0];
    std::array<Real, 2> eigenvalues;
    eigenvalues[0] = sum >= 0 ? .5*(sum + diff) : .5*(sum-diff);
    eigenvalues[1] = determinant/eigenvalues[0];
    if(eigenvalues[0] > eigenvalues[1]){
      std::swap(eigenvalues[0], eigenvalues[1]);
    }


    Real y = 2*S_[0][1], x = S_[0][0] - S_[1][1];
    Real theta = .5*atan2(y, x);


    auto R_theta = RotationMatrix_2x2(theta);
    auto R_theta_T = RotationMatrix_2x2(-theta);

    
    
    auto S_rotated_p = MultiplyMatrices_array(S_, R_theta);
    auto S_rotated = MultiplyMatrices_array(R_theta_T, S_rotated_p);

    
    
    array<Real, 2> eigenvalues_diagonalization{{S_rotated[0][0], S_rotated[1][1]}};
    if(eigenvalues_diagonalization[0] > eigenvalues_diagonalization[1]){
      std::swap(eigenvalues_diagonalization[0], eigenvalues_diagonalization[1]);
    }

    
    
    for(int j: {0, 1}){
      IsOk = IsOk &&
	abs(S_rotated[j][1-j]) < 10*numeric_limits<Real>::epsilon();
      IsOk = IsOk && abs(eigenvalues_diagonalization[j]-
			 eigenvalues[j]) <
	10*numeric_limits<Real>::epsilon();
    }

    
    
    return tuple{IsOk, S_rotated};
    
  }


}//namespace
