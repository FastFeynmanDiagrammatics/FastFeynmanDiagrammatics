namespace ffd::relation_vector{

  template<class T>
  auto
  WhichSet(T x,
	   Container<T> const& C){
    auto it = std::find(begin(C.container),
			end(C.container),
			x);
    if(it == end(C.container)) return std::array<uint, 2>{uint(-1), 0};
    std::size_t const dist = std::distance(begin(C.container), it);
    auto it2 = std::lower_bound(begin(C.split),
				end(C.split),
				dist);
    uint d0 = std::distance(begin(C.split), it2)
      - ((*it2) != dist);
    uint d1 = dist-C.split[d0];
    return std::array<uint, 2>{d0, d1};
  }

}//namespace
