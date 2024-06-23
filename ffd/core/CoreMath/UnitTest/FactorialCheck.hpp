

namespace ffd::core_math::unit_test{

  void FactorialCheck(){
    assert( Factorial(0) == 1 );
    assert( Factorial(1) == 1 );
    assert( Factorial(2) == 2 );
    assert( Factorial(3) == 6 );
    assert( Factorial(10) == 3628800 );
  }

}
