namespace ffd::class_tuple{

  template<template<typename> class, int, typename...>
  class ClassTupleInt;

  template<int j, typename T>
  struct ShellStruct{ T Content;};

  template<template<typename> class C, int j>
  class ClassTupleInt<C, j>{};

  template<template<typename> class C, int j, typename T, typename... Rest>
  class ClassTupleInt<C, j, T, Rest...>:
    public ShellStruct<j, C<T>>,
    public ClassTupleInt<C, j+1, Rest...>{};


  template<template<typename> class C, typename... T>
  using ClassTuple = ClassTupleInt<C, 0, T...>;
  
}//namespace ffd::class_tuple

