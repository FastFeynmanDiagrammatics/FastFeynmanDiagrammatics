namespace ffd::sigmadet_hubbard::unit_test{

  template < std::size_t n, class field_d = Real >
  void
  CDet_old_vs_new() {
    //    std::cerr << "CDet_old_vs_new\n";
    array2d<field_d, n, n> M;
    for ( ulong j = 0; j < n; ++j ) {
      for ( ulong k = 0; k < n; ++k ) {
	M(j, k) = (2*ffd::user_space::Proba()-1.)*(j!=k);
      } // for k in range(0, n)
    } // for j in range(0, n)
    ffd::vector3d<Real> Sigma ( n, n, (1ul<<n) );
    auto rho_T = Sigma;
    

    auto const bitmap = ffd::sigmadet_hubbard::CDet_Bitmap<n>();
    ffd::sigmadet_hubbard::CDet_Xi_r(M, Sigma, bitmap);
    
    auto [Xi_OLD, Z_OLD] = ffd::sigmadet_hubbard::compute_Xi_CDet_slow<n>(M);

    
    //    std::cerr << "Xi:" << std::endl;
    for(uint j=0; j<n; ++j){
      for(uint k=0; k<n; ++k){
        for(uint S=0; S<(1<<n); ++S){
	  //    std::cerr << j << " " << k << " " << S << "= " << Sigma(j, k, S)
	  //	          << " OLD= " << Xi_OLD[j+k*n][S] << " Xi" << std::endl;
    	  assert(( std::abs( Sigma(j, k, S) - Xi_OLD[j+k*n][S] ) < 1e-12 ));
        }
      }
    }

    
    ffd::sigmadet_hubbard::CDet_rho_T_r(M, Sigma, rho_T, bitmap);
    
    auto rho_OLD = ffd::sigmadet_hubbard::compute_rho_CDet_slow(M, Xi_OLD);

    
    //    std::cerr << "Xi:" << std::endl;
    for(uint j=0; j<n; ++j){
      for(uint k=0; k<n; ++k){
        for(uint S=0; S<(1<<n); ++S){
	  //    std::cerr << j << " " << k << " " << S << "= " << rho_T(k, j, S)
	  //	          << " OLD= " << rho_OLD[j+k*n][S] << " rho" << std::endl;
    	  assert(( std::abs( rho_T(k, j, S) - rho_OLD[j+k*n][S] ) < 1e-12 ));
        }
      }
    }

    
    ffd::sigmadet_hubbard::CDet_Sigma_r<n>(Sigma, rho_T, bitmap);
    auto Sigma_OLD = ffd::sigmadet_hubbard::compute_Sigma_CDet_slow(Xi_OLD, rho_OLD);


    //    std::cerr << "Sigma:" << std::endl;
    for(uint j=0; j<n; ++j){
      for(uint k=0; k<n; ++k){
        for(uint S=0; S<(1<<n); ++S){
	  //    std::cerr << j << " " << k << " " << S << "= " << Sigma(j, k, S)
	  //	          << " OLD= " << Sigma_OLD[j+k*n][S] << " Xi" << std::endl;
    	  assert(( std::abs( Sigma(j, k, S) - Sigma_OLD[j+k*n][S] ) < 1e-12 ));
        }
      }
    }

   
    
  }

}//namespace
