namespace ffd::four_point_vertex_hubbard{

  template<int order,
	   typename G0_t,
	   typename nilpoly_t>
  
  auto
  A_func(G0_t const& G0,
	 std::array<nilpoly_t, order*order*order*order> const& chi){
    std::array<nilpoly_t, order*order*order*order> a_func;
    a_func.fill(nilpoly_t(order));



    for(uint u1=0; u1<order; ++u1){
      for(uint d1=0; d1<order; ++d1){
	for(uint u2=0; u2<order; ++u2){
	  for(uint d2=0; d2<order; ++d2){
	    auto const prod = G0[u1*order+u2]*G0[d1*order+d2];
	    for(uint right=0; right<order*order; ++right){
	      a_func[(d1*order+u1)*order*order+right] +=
		prod*chi[(d2*order+u2)*order*order+right];
	    }
	  }
	}
      }
    }
    
    
    return a_func;
  }

}//namespace
