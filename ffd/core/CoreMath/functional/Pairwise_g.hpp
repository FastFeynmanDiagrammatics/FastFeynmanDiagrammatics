namespace ffd::functional{
  
  template<class T>
  auto&&
  Firstify(T&& x){
    return x.first;
  }

  template<class T>
  auto&&
  Secondify(T&& x){
    return x.second;
  }


  template<class func1_f,
	   class func2_f>
  auto
  Pairwise_g(func1_f const& f1,
	     func2_f const& f2){
    auto ret =
      [f1, f2]
      (auto&& ...x){
	return std::make_pair(f1(Firstify(std::forward<decltype(x)>(x))...),
			      f2(Secondify(std::forward<decltype(x)>(x))...));
      };
    return ret;
  }


  
}
