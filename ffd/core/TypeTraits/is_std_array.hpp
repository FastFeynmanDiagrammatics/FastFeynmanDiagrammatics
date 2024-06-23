namespace ffd::type_traits{

  template <typename T>
  struct is_std_array {
    static const bool value = false;
  };


  template <typename T, std::size_t S>
  struct is_std_array<std::array<T,S> > {
    static const bool value = true;
  };


  template<typename T>
  inline constexpr bool is_std_array_v = is_std_array<T>::value;

  

}//namespace
