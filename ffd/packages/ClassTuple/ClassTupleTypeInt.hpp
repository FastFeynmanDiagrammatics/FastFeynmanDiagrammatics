namespace ffd::class_tuple{

  template<template<typename, int> class, int, int, typename...>
  class ClassTupleTypeIntInt;

  template<template<typename, int> class C, int d, int j>
  class ClassTupleTypeIntInt<C, d, j>{};

  template<template<typename, int> class C, int d, int j, typename T, typename... Rest>
  class ClassTupleTypeIntInt<C, d, j, T, Rest...>:
    public ShellStruct<j, C<T, d>>,
    public ClassTupleTypeIntInt<C, d, j+1, Rest...>{
	};


  template<template<typename, int> class C, int d, typename... T>
  using ClassTupleTypeInt = ClassTupleTypeIntInt<C, d, 0, T...>;

}//namespace ffd::class_tuple
