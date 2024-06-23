

namespace ffd::eigenvalues_vectors{

  auto
  Diagonalize_2x2(std::array<std::array<Real, 2>, 2> const& S_){
    Real theta = 0.;

    
    assert((
    	    std::abs(S_[0][1] - S_[1][0]) <
    	    10*std::numeric_limits<Real>::epsilon()
    	    ));

    
    
    Real y = S_[0][1]+S_[1][0], x = -S_[0][0] + S_[1][1];
    return theta = .5*std::atan2(y, x);
  }


  

  auto
  Diagonalize_2x2(std::array<std::array<Complex, 2>, 2> const& H_){
    std::array<Real, 2> theta_phi{0., 0.};


    using ffd::core_math::Pi;
    
    
    auto& [theta, phi] = theta_phi;


    auto arg_H = -std::arg(H_[0][1]);

    phi = std::abs(arg_H) > .5*Pi ?
      arg_H - (arg_H>0)*Pi +(arg_H<0)*Pi :
      arg_H;



    std::array<std::array<Real, 2>, 2> S;
    for(int j: {0, 1}){
      S[j][j]   = std::real(H_[j][j]);
      S[j][1-j] = std::real( H_[j][1-j]*ffd::user_space::expI(phi*(1-2*j)) );
    }
    theta = Diagonalize_2x2(S);
    
    return theta_phi;
  }
  


}//namespace
