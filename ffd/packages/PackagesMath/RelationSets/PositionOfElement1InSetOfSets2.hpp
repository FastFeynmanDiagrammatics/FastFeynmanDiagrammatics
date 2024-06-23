namespace ffd::relation_sets{

  template<typename T>
  std::optional<int>
  PositionOfElement1InSetOfSets2(T const& element_,
				 std::set<std::set<T>> const& set_of_sets_){
    for(auto it = std::begin(set_of_sets_); it != std::end(set_of_sets_); ++it){
      std::optional<int> position = PositionOfElement1InSet2(element_, *it);
      if( position.has_value() ){
	return std::distance(std::begin(set_of_sets_), it);
      }
    }
    return {};
  }

}//namespace ffd::relation_sets
