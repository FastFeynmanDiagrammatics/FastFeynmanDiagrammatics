namespace ffd::functional{

  template<class func1_f,
	   class func2_f>
  auto
  PairFunctional_g(func1_f const& f1,
		   func2_f const& f2){
    auto ret =
      [f1, f2]
      (auto&& ...x){
	return std::make_pair(f1(std::forward<decltype(x)>...),
			      f1(std::forward<decltype(x)>...));
      };

    return ret;
  }
  
}
