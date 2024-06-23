

namespace ffd::integral_function::unit_test{

  void UnitTest(){
    using namespace std;
    IntegralFunction<Real> F([](Real x)->Real{return 1;}, 0);

    Real x = .2321;
    assert( abs(F(x)-x) < numeric_limits::epsilon() );

  }

}//namespace
