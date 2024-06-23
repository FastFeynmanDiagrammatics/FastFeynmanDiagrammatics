namespace ffd::four_point_vertex_hubbard{

  template<int order,
	   typename G0_t,
	   typename nilpoly_t>
  
  std::array<nilpoly_t, order*order*order*order>
  irreducible(G0_t const& G0,
	      std::array<nilpoly_t, order*order*order*order> const& a_func){
    std::array<nilpoly_t, order*order*order*order> T;
    T.fill(nilpoly_t(order));


    for( BinaryInt V = 1; V < two_n; ++V ){
      for(uint left=0; left<order*order; ++left){
	for(uint u1=0; u1<order; ++u1){
	  for(uint d1=0; d1<order; ++d1){
	    for(uint u2=0; u2<order; ++u2){
	      for(uint d2=0; d2<order; ++d2){
		for(BinaryInt S = ((V-1)&V); S != 0; S = ((S-1)&V)){
		  T[left*order*order+d1*order+u1][V] -=
		    T[left*order*order+d1*order+u1]
		Sigma[j*n+k][V] -= Sigma[j*n+l][S]*rho[l*n+k][V-S];
	      }

	     
	    }
	  }
	}
      }
    }
    return T;
  }


}//namespace
