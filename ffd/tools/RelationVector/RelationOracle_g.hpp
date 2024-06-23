namespace ffd::relation_vector{

  template<class T>
  auto 
  RelationOracle_g(std::vector<std::array<T, 2>> const& r){
    Container<T> C;
    for(auto x: r){
      AddRelation_r(x, C);
    }
    auto ret =
      [C]
      (auto... y)
      {
       if constexpr(sizeof...(y) == 0){
	   return C;
	 }else{
	 return WhichSet(y..., C);
       }
      };
    return ret;
  }

  
}//namespace
