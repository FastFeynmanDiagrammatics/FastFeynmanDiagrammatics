namespace ffd::functional{

  template<int j,
	   typename func_f>
  auto
  Ret_g(func_f const& f){
    auto ret =
      [f](auto&&... x){
	return std::get<j>(f(std::forward<decltype(x)>(x)...));
      };
    return ret;
  }

}
