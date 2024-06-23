

namespace ffd::core_math::unit_test{

  void assert_sqrt_int_test(int square, int result){
    assert( sqrt_int(square) == result );
  }

  
  void sqrt_int_test(){

    assert_sqrt_int_test(0, 0);

    assert_sqrt_int_test(1, 1);

    assert_sqrt_int_test(4, 2);


    assert_sqrt_int_test(9, 3);


    assert_sqrt_int_test(16, 4);

  }


}//namespace
