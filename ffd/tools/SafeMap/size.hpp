namespace ffd::user_space{

  template<typename value_t,
	   auto... keys>
  inline constexpr auto
  size(SafeMap<value_t, keys...> const&){
    return sizeof...(keys);
  }


}//namespace
