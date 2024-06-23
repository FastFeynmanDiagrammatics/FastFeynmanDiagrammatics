namespace ffd::sigmadet_hubbard{
  
  template<typename Xi_t=int,
	   typename rho_t=int>
  
  auto
  compute_Sigma_CDet_slow(Xi_t const& Xi,
			  rho_t const& rho){
    int const n = ffd::core_math::sqrt_int( size(Xi) );
    BinaryInt const two_n = (1<<n);
    auto Sigma = Xi;
    
    
    for( BinaryInt V = 1; V < two_n; ++V ){
      for( int j = 0; j < n; ++j ){
	auto const j_n = j*n;
	for( int k = 0; k < n; ++k ){
	  if( j != k ){
	    for( int l = 0; l < n; ++l ){
	      for( BinaryInt S = ((V-1)&V);
		   S != 0; S = ((S-1)&V) ){
		Sigma[j_n+k][V] -=
		  Sigma[j_n+l][S]*rho[l*n+k][V-S];
	      }
	    }
	  }
	}
      }
    }
    // for( BinaryInt V = 1; V < two_n; ++V ){
    //   for( int j = 0; j < n; ++j ){
    // 	for( int k = 0; k < n; ++k ){
    // 	  if( j != k ){
    // 	    for( int l = 0; l < n; ++l ){
    // 	      for( BinaryInt S = ((V-1)&V); S != 0; S = ((S-1)&V) ){
    // 		Sigma[j*n+k][V] -= Sigma[j*n+l][S]*rho[l*n+k][V-S];
    // 	      }
    // 	    }
    // 	  }
    // 	}
    //   }
    // }
    

    return Sigma;
  }

}//namespace
