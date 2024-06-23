

namespace ffd::packages_math::unit_test{

  struct test1{};


  template<typename T>
  struct test2{
    T y;


    constexpr auto how_many(){
      return sizeof_template_arguments(y);
    }

    
  };




  template<typename T1, typename T2>
  struct test3{
    T1 x;
    T2 y;
  };


  void test_sizeof_template_arguments(){

    test1 x;
    
    
    static_assert( sizeof_template_arguments(x) == 0);


    
    test2<test1> Y{};
    static_assert( sizeof_template_arguments(Y) == 1);

    
    static_assert( sizeof_template_arguments(Y.y) == 0);


    

    test2<test2<test1>> Z{};
    static_assert( sizeof_template_arguments(Z) == 1);


    static_assert( sizeof_template_arguments(Z.y) == 1);


    static_assert( sizeof_template_arguments(Z.y.y) == 0);




    
    test2<test3<test1, test2<test1>>> W{};

    static_assert( sizeof_template_arguments(W) == 1);
    
    static_assert( sizeof_template_arguments(W.y) == 2);

    static_assert( W.how_many() == 2);
    
  }


}//namespace
