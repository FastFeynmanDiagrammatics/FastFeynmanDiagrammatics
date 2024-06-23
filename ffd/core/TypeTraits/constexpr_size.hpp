namespace ffd::type_traits{


  template<typename T>
  struct constexpr_size{
    static const std::size_t value = 0;
  };
  

  
  template<typename T, std::size_t j>
  struct constexpr_size<std::array<T, j>>{
    static const std::size_t value = j;
  };


  template<typename T>
  inline constexpr std::size_t constexpr_size_v = constexpr_size<T>::value;

}//namespace
