

namespace ffd{

  template<template<typename> class C, typename... T>
  constexpr int size(ffd::class_tuple::ClassTuple<C, T...> const&){
    return sizeof...(T);
  }

  
  template<template<typename, int> class C, int d, typename... T>
  constexpr int size(ffd::class_tuple::ClassTupleTypeInt<C, d, T...> const&){
    return sizeof...(T);
  }
  
  
}//namespace ffd
