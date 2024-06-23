

namespace ffd::class_tuple::unit_test{

  template<typename T, int d>
  struct struct_merge_int{
    std::array<T, d> array;
  };

  
  void merge_type_int(){
    using ffd::merge;
    using ffd::size;
    
    ClassTupleTypeInt<struct_merge_int, 4, int, double, float> X1;
    ClassTupleTypeInt<struct_merge_int, 4, int> X2;
    ClassTupleTypeInt<struct_merge_int, 4> X3;

    auto Y1 = merge(X1, X2);
    assert( size(Y1) == 4);
    static_assert(std::is_same<decltype(Y1),
		  ClassTupleTypeInt<struct_merge_int, 4, int, double, float, int>>::value);
    
    auto Y2 = merge(X2, X1);
    assert( size(Y2) == 4);
    static_assert(std::is_same<decltype(Y2),
		  ClassTupleTypeInt<struct_merge_int, 4, int, int, double, float>>::value);

    auto Y3 = merge(X3, X1);
    assert( size(Y3) == 3);
    static_assert(std::is_same<decltype(Y3),
		  ClassTupleTypeInt<struct_merge_int, 4, int, double, float>>::value);
    
    auto Y4 = merge(X1, X3);
    assert( size(Y4) == 3);
    static_assert(std::is_same<decltype(Y4),
		  ClassTupleTypeInt<struct_merge_int, 4, int, double, float>>::value);

  }
  
}//namespace
