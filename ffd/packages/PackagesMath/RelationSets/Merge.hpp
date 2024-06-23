namespace ffd::relation_sets{

  template<typename T>
  std::set<T> Merge(std::set<T> const& set1_,
		    std::set<T> const& set2_){
    std::set<T> ret = set1_;
    for(auto const& element: set2_){
      if(ret.count(element) == 0){
	ret.insert(element);
      }
    }
    return ret;
  }

}//namespace ffd::relation_sets
