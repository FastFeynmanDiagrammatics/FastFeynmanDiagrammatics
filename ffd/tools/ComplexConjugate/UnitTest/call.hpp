

namespace ffd::complex_conjugate::unit_test{

  void UnitTest(){
    Real x = 1.;
    auto y = ComplexConjugate(x);
    static_assert( std::is_same<Real, decltype(y)>::value );

    int xx = 2;
    auto yy = ComplexConjugate(xx);
    assert( yy == xx );

    std::complex<int> z{1, 1};
    auto w = ComplexConjugate(z);
    static_assert( std::is_same<std::complex<int>, decltype(w)>::value );
    assert( std::imag(w) == -1 );
    

    Complex zz{1., 1.};
    auto ww = ComplexConjugate(zz);
    static_assert( std::is_same<Complex, decltype(ww)>::value );
    
  }

}//namespace
