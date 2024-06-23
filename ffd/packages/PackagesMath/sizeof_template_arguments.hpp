

namespace ffd::packages_math{

  template<template<auto...> typename T, auto... args>
  constexpr auto sizeof_template_arguments(T<args...> const&){
    return sizeof...(args);
  }


  template<template<typename...> typename T, typename... args>
  constexpr auto sizeof_template_arguments(T<args...> const&){
    return sizeof...(args);
  }


  template<typename T>
  constexpr auto sizeof_template_arguments(T const&){
    return 0;
  }


}//namespace
