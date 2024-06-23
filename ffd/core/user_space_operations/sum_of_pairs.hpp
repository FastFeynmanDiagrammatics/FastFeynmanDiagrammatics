namespace ffd::user_space{

  template<class T1, class T2>
  std::pair<T1, T2>
  operator+(std::pair<T1, T2> const& t0,
	    std::pair<T1, T2> const& t1){
    auto sum = std::make_tuple(t0.first, t0.second) +
      std::make_tuple(t1.first, t1.second);
    return std::make_pair(std::get<0>(sum),
			  std::get<1>(sum));
  }


}//namespace
