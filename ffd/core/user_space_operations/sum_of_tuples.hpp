namespace ffd::user_space{

  namespace operations{
    template<std::size_t r,
	     typename... T>
    void
    sum_of_tuples(std::tuple<T...> const& t0,
		  std::tuple<T...> const& t1,
		  std::tuple<T...>& tR){
      std::get<r>(tR) = std::get<r>(t0) + std::get<r>(t1);
      if constexpr(r+1 <sizeof...(T)) sum_of_tuples<r+1>(t0, t1, tR);
    }
  }

  template<typename... T>
  std::tuple<T...>
  operator+(std::tuple<T...> const& t0,
	    std::tuple<T...> const& t1){
    std::tuple<T...> tR;
    operations::sum_of_tuples<0>(t0, t1, tR);
    return tR;
  }

}//namespace
