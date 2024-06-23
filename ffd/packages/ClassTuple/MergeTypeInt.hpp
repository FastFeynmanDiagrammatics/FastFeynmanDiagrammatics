

namespace ffd{

  namespace class_tuple{
    
    using ffd::packages_math::IntType;
    
    template<IntType offset,
	     IntType first,
	     IntType... last,
	     template<typename, int> class C,
	     int d,
	     typename... T1,
	     typename... T2>
    void PartiallyCopy1In2(ClassTupleTypeInt<C, d, T1...> const& X1,
			   ClassTupleTypeInt<C, d, T2...>& X2){
      using ffd::get;
      ( ( get<first+offset>(X2) = get<first>(X1) ) ,
	... , ( get<last+offset>(X2) = get<last>(X1) ) );
    }


    template<IntType offset,
	     IntType... last,
	     template<typename, int> class C,
	     int d,
	     typename... T1,
	     typename... T2>
    void PartiallyCopy1In2(ClassTupleTypeInt<C, d, T1...> const& X1,
			   ClassTupleTypeInt<C, d, T2...>& X2,
			   ffd::packages_math::StaticRangeSequence<0, 0, last...>){
      PartiallyCopy1In2<offset, last...>(X1, X2);
    }

  }

  
  template<template<typename, int> class C,
	   int d,
	   typename... T1,
	   typename... T2>
  auto merge(ffd::class_tuple::ClassTupleTypeInt<C, d, T1...> const& X1,
	     ffd::class_tuple::ClassTupleTypeInt<C, d, T2...> const& X2){
    ffd::class_tuple::ClassTupleTypeInt<C, d, T1..., T2...> ret;
    if constexpr(sizeof...(T1) > 0){
	ffd::packages_math::StaticRangeSequence<0, sizeof...(T1)> S1;
	ffd::class_tuple::PartiallyCopy1In2<0>(X1, ret, S1);
      }
    if constexpr(sizeof...(T2) > 0){
	ffd::packages_math::StaticRangeSequence<0, sizeof...(T2)> S2;
	ffd::class_tuple::PartiallyCopy1In2<sizeof...(T1)>(X2, ret, S2);
      }
    return ret;
  }


}//namespace ffd::class_tuple
