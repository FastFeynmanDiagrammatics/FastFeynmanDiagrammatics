namespace ffd::math_tools::unit_test{

  void UnitTest(){

    static_assert(( log2_int(1) == 0 ));

    static_assert(( log2_int( 3 )  ==  1 ));

    static_assert(( log2_int( 4 )  ==  2 ));

    static_assert(( log2_int( 20 )  ==  4 ));

    mean_diff();

//    test_print_matrix();
  }

}//namespace
