namespace ffd::type_traits{

  template<typename T>
  struct element_of_container{
    using type = typename std::decay_t<decltype(std::declval<T>()[0])>;
  };


  template<typename T>
  using element_of_container_t = typename element_of_container<T>::type;


}//namespace
