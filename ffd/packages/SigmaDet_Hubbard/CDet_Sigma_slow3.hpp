namespace ffd::sigmadet_hubbard{

  template<std::size_t n, class field_d = Real>
  void
  CDet_Sigma_slow3_r(array3d<field_d, n, n, (1ul<<n)>& __restrict__ Xi,
		     array3d<field_d, n, n, (1ul<<n)> const&  __restrict__ rho_T
		     ){
    static_assert(( n <= 14 ));
    ulong constexpr two_n = (1<<n);
    for( ulong V = 1; V < two_n; ++V ){
      for( ulong S = ((V-1)&V); S != 0; S = ((S-1)&V)){
        for( ulong j = 0; j < n; ++j ){
          if( (S&(1ul<<j)) != 0){
            auto const S_n_j_n = n*n*S+j*n;
            ulong k=0;
            for(; k < j; ++k ){
              if( ((V-S)&(1ul<<k)) != 0){
                auto const V_S_n_k = n*n*(V-S)+n*k;
		field_d acc = 0;
#pragma unroll
		for( ulong l = 0; l < bm(0, S); ++l ){
                  acc += Xi[S_n_j_n+bm(l, S)] * rho_T[V_S_n_k+bm(l, S)];
                }
		Xi[V_n_j_n_k] = acc;
		if constexpr ( n > 1 ) {
		    acc += Xi[S_n_j_n+1] * rho_T[V_S_n_k+1];
		  }
		if constexpr ( n > 2 ) {
		    acc += Xi[S_n_j_n+2] * rho_T[V_S_n_k+2];
		  }
		if constexpr ( n > 3 ) {
		    acc += Xi[S_n_j_n+3] * rho_T[V_S_n_k+3];
		  }
		if constexpr ( n > 4 ) {
		    acc += Xi[S_n_j_n+4] * rho_T[V_S_n_k+4];
		  }
		if constexpr ( n > 5 ) {
		    acc += Xi[S_n_j_n+5] * rho_T[V_S_n_k+5];
		  }
		if constexpr ( n > 6 ) {
		    acc += Xi[S_n_j_n+6] * rho_T[V_S_n_k+6];
		  }
		if constexpr ( n > 7 ) {
		    acc += Xi[S_n_j_n+7] * rho_T[V_S_n_k+7];
		  }
		if constexpr ( n > 8 ) {
		    acc += Xi[S_n_j_n+8] * rho_T[V_S_n_k+8];
		  }
		if constexpr ( n > 9 ) {
		    acc += Xi[S_n_j_n+9] * rho_T[V_S_n_k+9];
		  }
		if constexpr ( n > 10 ) {
		    acc += Xi[S_n_j_n+10] * rho_T[V_S_n_k+10];
		  }
		if constexpr ( n > 10 ) {
		    acc += Xi[S_n_j_n+10] * rho_T[V_S_n_k+10];
		  }
		if constexpr ( n > 11 ) {
		    acc += Xi[S_n_j_n+11] * rho_T[V_S_n_k+11];
		  }
		if constexpr ( n > 12 ) {
		    acc += Xi[S_n_j_n+12] * rho_T[V_S_n_k+12];
		  }
		if constexpr ( n > 13 ) {
		    acc += Xi[S_n_j_n+13] * rho_T[V_S_n_k+13];
		  }
		
		Xi[n*(n*V+j)+k] -= acc;
              }
            }
            ++k;
            for(; k < n; ++k ){
              if( ((V-S)&(1ul<<k)) != 0){
                auto const V_S_n_k = n*n*(V-S)+n*k;
		field_d acc = Xi[S_n_j_n] * rho_T[V_S_n_k];
// #pragma unroll
//                 for( ulong l = 0; l < n; ++l ){
//                   acc += Xi[S_n_j_n+l] * rho_T[V_S_n_k+l];
//                 }
		if constexpr ( n > 1 ) {
		    acc += Xi[S_n_j_n+1] * rho_T[V_S_n_k+1];
		  }
		if constexpr ( n > 2 ) {
		    acc += Xi[S_n_j_n+2] * rho_T[V_S_n_k+2];
		  }
		if constexpr ( n > 3 ) {
		    acc += Xi[S_n_j_n+3] * rho_T[V_S_n_k+3];
		  }
		if constexpr ( n > 4 ) {
		    acc += Xi[S_n_j_n+4] * rho_T[V_S_n_k+4];
		  }
		if constexpr ( n > 5 ) {
		    acc += Xi[S_n_j_n+5] * rho_T[V_S_n_k+5];
		  }
		if constexpr ( n > 6 ) {
		    acc += Xi[S_n_j_n+6] * rho_T[V_S_n_k+6];
		  }
		if constexpr ( n > 7 ) {
		    acc += Xi[S_n_j_n+7] * rho_T[V_S_n_k+7];
		  }
		if constexpr ( n > 8 ) {
		    acc += Xi[S_n_j_n+8] * rho_T[V_S_n_k+8];
		  }
		if constexpr ( n > 9 ) {
		    acc += Xi[S_n_j_n+9] * rho_T[V_S_n_k+9];
		  }
		if constexpr ( n > 10 ) {
		    acc += Xi[S_n_j_n+10] * rho_T[V_S_n_k+10];
		  }
		if constexpr ( n > 10 ) {
		    acc += Xi[S_n_j_n+10] * rho_T[V_S_n_k+10];
		  }
		if constexpr ( n > 11 ) {
		    acc += Xi[S_n_j_n+11] * rho_T[V_S_n_k+11];
		  }
		if constexpr ( n > 12 ) {
		    acc += Xi[S_n_j_n+12] * rho_T[V_S_n_k+12];
		  }
		if constexpr ( n > 13 ) {
		    acc += Xi[S_n_j_n+13] * rho_T[V_S_n_k+13];
		  }
		Xi[n*(n*V+j)+k] -= acc;
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
