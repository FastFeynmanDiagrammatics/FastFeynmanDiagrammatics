namespace ffd::sigmadet_hubbard{

  template<std::size_t n, class field_d = Real>
  void
  CDet_rho_slow_r(
		  array2d<field_d, n, n> const& __restrict__ G0,
		  array3d<field_d, n, n, (1ul<<n)>& __restrict__ Xi,
		  array3d<field_d, n, n, (1ul<<n)>& __restrict__ rho
		  ){
    rho.fill(0);
    for ( ulong S=1; S<(1ul<<n); ++S ) {
      for ( ulong j = 0; j < n; ++j ) {
	for ( ulong k = 0; k < n; ++k ) {
	  if ( j != k ) {
	    for ( ulong l = 0; l < n; ++l ) {
	      rho(k, j, S)
		 +=
		G0(l, j)*
		Xi(k, l, S);
	    } // for l in range(0, n)
	  } // if j != k 
	} // for k in range(0, n)
      } // for j in range(0, n)
    } // for S in range(1, (1ul<<n))
    return;
  } // function CDet_rho_r

} // namespace
