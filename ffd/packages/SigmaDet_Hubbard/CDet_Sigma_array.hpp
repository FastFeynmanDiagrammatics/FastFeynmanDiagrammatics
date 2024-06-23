namespace ffd::sigmadet_hubbard{

  template<std::size_t n, class field_d = Real>
  void
  CDet_Sigma_array_r(
		     array3d<field_d, n, n, (1ul<<n)>& __restrict__ Sigma,
		     array3d<field_d, n, n, (1ul<<n)> const&  __restrict__ rho_T,
		     array2d<std::uint8_t, n+1, (1ul<<n)> const& __restrict__ bm
		     ){
    ulong constexpr two_n = (1<<n);
    ulong constexpr n2 = n*n;
    for ( ulong V = 1; V < two_n; ++V ) {
      ulong const V_n = n2*V;
      for ( ulong S = ((V-1)&V); S != 0; S = ((S-1)&V)) {
	ulong const bm_0_S = bm(0, S);
	ulong const bm_0_V_S = bm(0, V-S);
	ulong const V_S_n = n2*(V-S);
	ulong const S_n = n2*S;
	for ( ulong j = 0; j < bm_0_S; ++j ) {
	  ulong const b_n_j = n*bm(j+1, S);
	  ulong const S_n_j_n = S_n + b_n_j;
	  ulong const V_n_j = V_n + b_n_j;
	  for ( ulong k = 0; k < bm_0_V_S; ++k ) {
	    ulong const b_k = bm(k+1, V-S);
	    ulong const V_S_n_k = V_S_n + n*b_k;
	    field_d acc = 0;
#pragma GCC unroll 4
#pragma GCC ivdep
	    for ( ulong l = 0; l < bm_0_S; ++l ) {
	      ulong const b_l = bm(l+1, S);
	      acc += Sigma[S_n_j_n+b_l] *
		rho_T[V_S_n_k+b_l];
	    } // for l in range(0, bm(0, S))
	    Sigma[V_n_j+b_k] -= acc;
	  }
	}
      }
    }
    return;
  } // routine CDet_Sigma

}//namespace
