

namespace ffd::class_tuple::unit_test{

  template<typename T>
  struct struct_test{
    T element = 1;
  };

  template<typename T>
  void do_nothing(T){}

  void class_tuple(){
    ClassTuple<struct_test, int, double, char> X;
    do_nothing(X);
  }

}//namespace
