namespace ffd::normal_principal_blocks{
  
  template<typename coord_t, typename field_t>
  std::vector<std::pair<std::vector<QuantumFieldPair<coord_t>>,
			std::array<std::vector<int>, 2>>>
  Blocks(std::vector<QuantumFieldSum<coord_t, field_t>> const& I,
	 QuantumFieldSum<coord_t, field_t> const& V_ext,
	 ffd::relation_sets::RelationSets<QF> const& corr_blocks){
    std::vector<std::pair<std::vector<QuantumFieldPair<coord_t>>,
			  std::array<std::vector<int>, 2>>>
      ret(size(corr_blocks));
    std::size_t const order = size(I);
    std::size_t const n_ext = size(V_ext.fields[0])/2;
    assert(((size(V_ext.fields[0])/2)*2 == size(V_ext.fields[0])));
    auto I_S = I[0];
    for(uint j=1; j<order; ++j){
      I_S *= I[j];
    }
    I_S *= V_ext;
    
    
    for(uint j=0; j<order; ++j){
      auto p0 = fill_g0_matrix::pop_front_qfpair(I_S);
      auto p1 = fill_g0_matrix::pop_front_qfpair(I_S);
      // std::cerr << p0.fields[0].first << 'x' <<
      // 	p0.fields[1].first << '\n';
      // std::cerr << p1.fields[0].first << 'x' <<
      // 	p1.fields[1].first << '\n';
      if(!ffd::relation_sets::
	 Are12InRelation3(p0.fields[0].first,
			  Bar(p0.fields[1].first),
			  corr_blocks).value() ||
	 !ffd::relation_sets::
	 Are12InRelation3(p1.fields[0].first,
			  Bar(p1.fields[1].first),
			  corr_blocks).value()){
	std::swap(p0.fields[0], p1.fields[0]);
	if(p0.fields[0].first.statistics ==
	   phys::fermi){
	  p0.sign = -p0.sign;
	}
      }
      assert((ffd::relation_sets::
	      Are12InRelation3(p0.fields[0].first,
			       Bar(p0.fields[1].first),
			       corr_blocks).value()));
      assert((ffd::relation_sets::
	      Are12InRelation3(p1.fields[0].first,
			       Bar(p1.fields[1].first),
			       corr_blocks).value()));
      
      
      auto block0 = ffd::relation_sets::
	PositionOfElement1InSetOfSets2(p0.fields[0].first, corr_blocks).value();
      auto block1 = ffd::relation_sets::
	PositionOfElement1InSetOfSets2(p1.fields[0].first, corr_blocks).value();
      
      
      ret[block0].first.push_back(p0);
      ret[block0].second[0].push_back(j);
      ret[block1].first.push_back(p1);
      ret[block1].second[0].push_back(j);
    }

    
    for(int j=0; j<n_ext; ++j){
      auto p0 = fill_g0_matrix::pop_front_qfpair(I_S);
      auto block0 = ffd::relation_sets::
	PositionOfElement1InSetOfSets2(p0.fields[0].first, corr_blocks).value();
      ret[block0].first.push_back(p0);
      ret[block0].second[1].push_back(j);
    }

    
    return ret;
  }

}//namespace
