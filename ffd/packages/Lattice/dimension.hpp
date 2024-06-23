namespace ffd::lattice{

  
  template<template<int, auto...> typename T, int d, auto... args>
  constexpr int dimension(T<d, args...> const&){
    return d;
  }
  

}//namespace
