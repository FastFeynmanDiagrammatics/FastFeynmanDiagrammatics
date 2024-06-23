namespace ffd::user_space::fill_g0_matrix{
  
  template<typename coord_t, typename field_t>
  QuantumFieldPair<coord_t>
  pop_front_qfpair(QuantumFieldSum<coord_t, field_t>& qfs){
    int sign = 1;
    std::pair<QF, coord_t> psi;
    auto& qf_prod = qfs.fields[0];
    for(uint j=0; j<size(qf_prod); ++j){
      if(qf_prod[j].first.direction ==
	 phys::in){
	if(j%2 == 0 &&
	   qf_prod[j].first.statistics ==
	   phys::fermi){
	  sign = - sign;
	}
	psi = qf_prod[j];
	for(uint u=j; u<size(qf_prod)-1; ++u){
	  qf_prod[u] = qf_prod[u+1];
	}
	break;
      }
    }
    qf_prod.resize(size(qf_prod)-1);
    
    
    std::pair<QF, coord_t> bpsi;
    for(uint j=0; j<size(qf_prod); ++j){
      if(qf_prod[j].first.direction ==
	 phys::ou){
	if(j%2 == 0 &&
	   qf_prod[j].first.statistics ==
	   phys::fermi){
	  sign = - sign;
	}
	bpsi = qf_prod[j];
	for(uint u=j; u<size(qf_prod)-1; ++u){
	  qf_prod[u] = qf_prod[u+1];
	}
	break;
      }
    }
    qf_prod.resize(size(qf_prod)-1);
    
	
    auto ret = QuantumFieldPair<coord_t>(psi, bpsi);
    ret.sign *= sign;
    return ret;
  }

}//namespace
