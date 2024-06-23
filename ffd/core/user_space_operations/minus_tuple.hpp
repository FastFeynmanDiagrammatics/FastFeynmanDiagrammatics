namespace ffd::user_space{

  namespace operations{
    template<std::size_t r,
	     typename... T>
    void
    minus_tuple(
		std::tuple<T...>& tR){
      std::get<r>(tR) =  - std::get<r>(tR);
      if constexpr(r+1 < sizeof...(T)) minus_tuple<r+1>(tR);
    }
  }

  
  template<typename... T>
  std::tuple<T...>
  operator-(std::tuple<T...> const& t0
	    ){
    std::tuple<T...> tR  = t0;
    operations::minus_tuple<0>(tR);
    return tR;
  }

}//namespace
