namespace ffd{
  
  template<auto element,
	   typename value_t,
	   std::size_t size_v>
  value_t&
  get(ffd::user_space::SafeArray<value_t, size_v>& Array_){
    static_assert( std::size_t(element) >= 0 , "SafeArray out of bounds");
    static_assert( std::size_t(element) < size_v, "SafeArray out of bounds");
    
    
    return Array_.Array[element];
  }



  template<auto element,
	   typename value_t,
	   std::size_t size_v>
  value_t const&
  get(ffd::user_space::SafeArray<value_t, size_v> const& Array_){
    static_assert( std::size_t(element) >= 0 , "SafeArray out of bounds");
    static_assert( std::size_t(element) < size_v, "SafeArray out of bounds");
    
    
    return Array_.Array[element];
  }

  
}//namespace
