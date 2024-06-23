namespace ffd::sigmadet_hubbard{

  template<std::size_t n, class field_d = Real>
  void
  CDet_Sigma_slow_r(array3d<field_d, n, n, (1ul<<n)>& __restrict__ Xi,
		    array3d<field_d, n, n, (1ul<<n)> const&  __restrict__ rho){
    ulong constexpr two_n = (1<<n);
    ulong constexpr nn = n*n;
    for( ulong V = 1; V < two_n; ++V ){
      auto const V_n = nn*V;
      for( ulong S = ((V-1)&V); S != 0; S = ((S-1)&V)){
        auto const V_S = V-S;
        auto const V_S_n = nn*V_S;
        auto const S_n = nn*S;
        for( ulong j = 0; j < n; ++j ){
          if( (S&(1<<j)) != 0){
            auto const j_n = j*n;
            auto const V_n_j_n = V_n+j_n;
            auto const S_n_j_n = S_n+j_n;
            ulong k=0;
            for(; k < j; ++k ){
              if( (V_S&(1<<k)) != 0){
                auto const V_S_n_k = V_S_n+k;
                auto const V_n_j_n_k = V_n_j_n+k;
#pragma unroll
                for( ulong l = 0; l < n; ++l ){
                  Xi[V_n_j_n_k] -= Xi[S_n_j_n+l]*rho[V_S_n_k+l*n];
                }
              }
            }
            ++k;
            for(; k < n; ++k ){
              if( (V_S&(1<<k)) != 0){
                auto const V_S_n_k = V_S_n+k;
                auto const V_n_j_n_k = V_n_j_n+k;
#pragma unroll
                for( ulong l = 0; l < n; ++l ){
                  Xi[V_n_j_n_k] -= Xi[S_n_j_n+l]*rho[V_S_n_k+l*n];
                }
              }
            }
          }
        }
      }
    }

    // for( ulong V=1; V<two_n; ++V ){
    //   for( ulong S=((V-1)&V); S!=0; S=((S-1)&V) ){
    //     for( ulong j=0; j<n; ++j ){
    //       if( (S&(1<<j))!=0 ){
    //         for(ulong k=0; k<n; ++k ){
    //           if( k!= j ){
    //             if( (V-S&(1<<k)) != 0){
    //               for( ulong l = 0; l < n; ++l ){
    //                 Xi[n*n*V+j*n+k] -= Xi[n*n*S+j*n+l]*rho[V-S*n*n+l*n+k];
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    // }
    return;
  } // routine CDet_Sigma

}//namespace
