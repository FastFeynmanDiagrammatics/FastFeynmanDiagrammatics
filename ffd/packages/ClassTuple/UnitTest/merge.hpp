

namespace ffd::class_tuple::unit_test{
  template<typename T>
  struct struct_merge{
    T element;
  };

  void merge(){
    using ffd::merge;
    using ffd::size;
    
    ClassTuple<struct_merge, int, double, float> X1;
    ClassTuple<struct_merge, int> X2;
    ClassTuple<struct_merge> X3;

    auto Y1 = merge(X1, X2);
    assert( size(Y1) == 4);
    static_assert(std::is_same<decltype(Y1),
		  ClassTuple<struct_merge, int, double, float, int>>::value);
    
    auto Y2 = merge(X2, X1);
    assert( size(Y2) == 4);
    static_assert(std::is_same<decltype(Y2),
		  ClassTuple<struct_merge, int, int, double, float>>::value);

    auto Y3 = merge(X3, X1);
    assert( size(Y3) == 3);
    static_assert(std::is_same<decltype(Y3),
		  ClassTuple<struct_merge, int, double, float>>::value);
    
    auto Y4 = merge(X1, X3);
    assert( size(Y4) == 3);
    static_assert(std::is_same<decltype(Y4),
		  ClassTuple<struct_merge, int, double, float>>::value);
        
  }

}//namespace
