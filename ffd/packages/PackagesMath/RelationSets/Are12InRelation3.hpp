namespace ffd::relation_sets{

  template<typename T>
  std::optional<bool>
  Are12InRelation3(T const& e0,
		   T const& e1,
		   std::set<std::set<T>> const& rel_set){
    auto pos0 = PositionOfElement1InSetOfSets2(e0, rel_set);
    auto pos1 = PositionOfElement1InSetOfSets2(e1, rel_set);
    if(!pos0.has_value() ||
       !pos1.has_value()){
      return {};
    }
    return pos0.value() == pos1.value();
  }

}
