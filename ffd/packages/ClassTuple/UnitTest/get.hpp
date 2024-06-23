

namespace ffd::class_tuple::unit_test{

  template<typename T>
  struct struct_get{
    T element;
  };

  template<typename T, int j>
  struct struct_get_int{
    std::array<T, j> array;
  };

  
  void get(){
    using ffd::get;
    ClassTuple<struct_get, int, char, char> X1;

    get<0>(X1).element = 3;
    assert(get<0>(X1).element == 3);

    get<1>(X1).element = -2;
    assert(get<1>(X1).element == -2);
    
    get<2>(X1).element = 'a';
    assert(get<2>(X1).element == 'a');


    ClassTupleTypeInt<struct_get_int, 2, int, char, char> X2;

    std::array<int, 2> array0 = {1, 2};
    get<0>(X2).array = array0;
    assert(get<0>(X2).array == array0);

    std::array<char, 2> array2 = {'d', 'e'};
    get<2>(X2).array = array2;
    assert(get<2>(X2).array == array2);
    
  }

}//namespace
