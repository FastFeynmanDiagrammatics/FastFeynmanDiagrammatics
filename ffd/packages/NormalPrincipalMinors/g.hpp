namespace ffd::user_space{
    
  template<int order,
	   class nilpoly_d,
	   class coord_d, class field_d,
	   class G0_f,
	   class b_oracle_f,
	   class alpha_shift_f>
  std::array<nilpoly_d, 2>
  NormalPrincipalMinors(std::array<QFSum<coord_d, field_d>, order> const& I,
			QFSum<coord_d, field_d> const& V_ext,
			G0_f const& G0,
			b_oracle_f const& b_oracle,
			a_shift_f const& a_shift){
    std::array<nilpoly_d, 2> ret;
    ret.fill(nilpoly_d(int(order)));
    for(BinaryInt S=0; S<(1<<order); ++S) for(int s: {0, 1}) ret[s][S] = 1.;

    std::vector<std::array<std::pair<QField, coord_d>
    auto chis = create_chis()
    
    for(uint j=0; j<size(vec_blocks); ++j){
      auto chis = vec_blocks[j].first;
      auto int_vtx = vec_blocks[j].second[0];
      auto ext_vtx = vec_blocks[j].second[1];


      timer_v.ini();
      auto M = fill_G0_matrix(G0, chis);
      timer_v.fin();

      //      ffd::math_tools::print_matrix<3>(M);
      if constexpr(!std::is_same_v<int, alpha_shift_t>){
	  M = alpha_shift(M, order);
	}
      // ffd::math_tools::print_matrix<3>(M);
      
      
      int step_size = 1;
      if(contains_double_element(int_vtx)){
	step_size = 2;
	for(uint k=0; k<size(int_vtx); ++k){
	  if(!is_double_element(k, int_vtx)){
	    pad_matrix(M, 2*k+1);
	  }
	}
	if((size(ext_vtx)%2) == 1){
	  pad_matrix(M);
	}
	remove_double_elements(int_vtx);
      }


      //      timer_v.ini();
      // auto PM = ffd::principal_minors::
      // 	DeterminantPrincipalMinors(M, step_size);
      auto PM = ffd::principal_minors::
	PrincipalMinors(M);
      //timer_v.fin();
      
      // for(BinaryInt S=0; S<size(PM); ++S){
      // 	std::cerr << S << ' ' << PM[S] << '\n';
      // }

      
      BinaryInt v = (1<<size(int_vtx))-1;
      BinaryInt v_ext = size(PM) - 1 - v;
      BinaryInt V = 0;
      for(uint j=0; j<size(int_vtx); ++j){
	V += (1<<(int_vtx[j]));
      }
      for(BinaryInt S = V, s = v;
	  S > 0;
	  S = ((S-1)&V), --s){
	ret[0][S] *= PM[s];
	ret[1][S] *= PM[s+v_ext];
      }
      ret[0][0] *= PM[0];
      ret[1][0] *= PM[v_ext];
    }


    // std::cerr << ret[0] << '\n';
    // std::cerr << ret[1] << '\n';

    
    std::vector<field_t> interaction_strength;
    for(uint j=0; j<order; ++j){
      interaction_strength.push_back(-I[j].coef[0]);
    }
    for(int s: {0, 1}) ret[s] = ret[s](interaction_strength);
    
    
    return ret;
  }


}//namespace
