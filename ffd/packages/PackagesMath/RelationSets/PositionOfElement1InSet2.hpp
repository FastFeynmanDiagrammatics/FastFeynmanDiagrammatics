namespace ffd::relation_sets{

  template<typename T>
  std::optional<int>
  PositionOfElement1InSet2(T const& element_,
			   std::set<T> const& set_){
    auto position = set_.find(element_);
    if( position == std::end(set_) ){
      return {};
    }
    return std::distance(std::begin(set_), position);
  }
  
}//namespace ffd::relation_sets
