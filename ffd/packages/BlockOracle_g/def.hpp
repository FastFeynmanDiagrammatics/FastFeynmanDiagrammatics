namespace ffd::user_space{

  template<class coord_d, class field_d>
  auto
  BlockOracle_g(QFSum<coord_d, field_d> const& H0){
    assert((size(H0.fields)%2 == 0));
    std::vector<std::array<QField, 2>> relation_vec;
    for(uint term=0; term<size(H0.fields)/2; ++term){
      auto f0 = H0.fields[2*term].first;
      if(Direction(f0) == phys::ou) f0 = Bar(f0);
      auto f1 = H0.fields[2*term+1].first;
      if(Direction(f1) == phys::ou) f1 = Bar(f1);
      relation_vec.push_back({f0, f1});
    }
    auto oracle = relation_vector::
      RelationOracle_g(relation_vec);
    auto flip_if_need =
      [](QField q){
	if(Direction(q) == phys::ou) q = Bar(q);
	return q;
      };
    auto ret =
      [oracle, flip_if_need]
      (auto... q){
	if constexpr(sizeof...(q) == 0){
	    return oracle();
	  }else{
	  return oracle(flip_if_need(q...));
	}
      };
    return ret;
  }

}//namespace
