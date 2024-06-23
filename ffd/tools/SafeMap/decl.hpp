namespace ffd::user_space{
  
  template<typename value_t,
	   auto... keys>
  struct SafeMap{
    std::array<value_t, sizeof...(keys)> values;

    SafeMap() {}

    template<template<auto...> typename int_seq_t>
    SafeMap(int_seq_t<keys...>, std::array<value_t, sizeof...(keys)> values_):
      values(values_) {}

    
    template<typename key_t>
    value_t& operator[](key_t key){
      int position = ffd::user_space::safe_map::
	dynamic_position_of_value_in_sequence<key_t, keys...>(key);
      assert(( position != -1 ));
      return values[position];
    }


    template<typename key_t>
    const value_t& operator[](key_t key) const{
      int position = ffd::user_space::safe_map::
	dynamic_position_of_value_in_sequence<key_t, keys...>(key);
      assert(( position != -1 ));
      return values[position];
    }

  
  };

}//namespace


