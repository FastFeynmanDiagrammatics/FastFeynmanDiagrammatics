

namespace ffd::class_tuple::unit_test{

  template<typename T, int j>
  struct struct_test_int{
    T element = j;
  };

  
  void class_tuple_type_int(){
    ClassTupleTypeInt<struct_test_int, 3, int, double, char, char, int> X;
    do_nothing(X);
  }
  

}//namespace
