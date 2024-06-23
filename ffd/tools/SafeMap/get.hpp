namespace ffd{

  template<auto key, typename value_t, auto... keys>
  
  value_t const&
  get(ffd::user_space::SafeMap<value_t, keys...> const& M){
    constexpr int position = ffd::user_space::safe_map::
      static_position_of_value_in_sequence<key, keys...>();
    static_assert( position != -1 , "SafeMap: key not found");
    return M.values[position];
  }



  template<auto key, typename value_t, auto... keys>
  
  value_t&
  get(ffd::user_space::SafeMap<value_t, keys...>& M){
    return const_cast<value_t&>(ffd::get<key>(static_cast<ffd::user_space::SafeMap<value_t, keys...> const &>(M)));
  }

}//namespace 
