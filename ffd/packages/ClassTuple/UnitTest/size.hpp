

namespace ffd::class_tuple::unit_test{


  template<typename T>
  struct struct_test_size{
    T x1, x2;
  };

  
  template<typename T, int j>
  struct struct_test_size_int{
    T x1, x2;
    std::array<char, j> x3;
  };

  
  void size(){
    using ffd::size;
    ClassTuple<struct_test_size, double, int, char, char> X1;
    assert(size(X1) == 4);

    ClassTupleTypeInt<struct_test_size_int, 5, double, float,int, char, char> X2;
    assert(size(X2) == 5);

  }

}//namespace
