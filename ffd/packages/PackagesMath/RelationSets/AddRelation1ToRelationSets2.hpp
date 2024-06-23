namespace ffd::relation_sets{

  template<typename T>
  void
  AddRelation1ToRelationSets2(Element0IsInRelationWithElement1<T>
			      const& relation_elements_,
			      RelationSets<T>& relation_sets_){
    assert( std::size(relation_elements_) == 2 );
    std::array<std::optional<int>, 2> positions;
    for(auto j: {0, 1}){
      positions[j] = PositionOfElement1InSetOfSets2(relation_elements_[j], relation_sets_);
    }
    if(!positions[0].has_value() && !positions[1].has_value()){
      std::set<T> set_to_add{relation_elements_[0], relation_elements_[1]};
      relation_sets_.insert(set_to_add);
    }else{
      for(auto j: {0, 1}){
	if(positions[j].has_value() && !positions[1-j].has_value()){
	  auto it = std::begin(relation_sets_);
	  for(int k=0; k < positions[j].value(); ++k){
	    ++it;
	  }
	  auto set_j = *it;
	  set_j.insert(relation_elements_[1-j]);
	  relation_sets_.erase(it);
	  relation_sets_.insert(set_j);
	}
      }
      if(positions[0].has_value() && positions[1].has_value()){
	if(positions[0].value() != positions[1].value()){
	  std::array<std::set<T>, 2> sets_0_1;
	  for(auto j: {0, 1}){
	    positions[j] = PositionOfElement1InSetOfSets2(relation_elements_[j], relation_sets_);
	    auto it = std::begin(relation_sets_);
	    for(int k=0; k < positions[j].value(); ++k){
	      ++it;
	    }
	    sets_0_1[j] = *it;
	    relation_sets_.erase(it);
	  }
	  relation_sets_.insert(Merge(sets_0_1[0], sets_0_1[1]));
	}
      }
    }
  }

}//namespace ffd::relation_sets
