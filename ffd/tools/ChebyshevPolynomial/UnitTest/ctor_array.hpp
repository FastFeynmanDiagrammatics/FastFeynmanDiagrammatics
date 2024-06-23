namespace ffd::chebyshev_polynomial::unit_test{
  
  void ctor_array(){
    
    constexpr auto N_Cheby = 30;
    Real Beta = 5.;

    std::array<Real, N_Cheby> values;
    auto f = [](Real x){return x*x*sin(x);};
    for( int j=0; j<N_Cheby; ++j){
      Real tau = .5*Beta*(1+cos(ffd::core_math::Pi*(j+.5)/N_Cheby));
      values[j] = f( tau );
    }


    auto F = ChebyshevPolynomial<Real>(values, {0., Beta});

    assert(( std::abs(F(0)-f(0)) < 1e-10 ));
    // std::cerr<<F(0)<<" "<<f(0)<<std::endl;
    // std::cerr<<F(Beta)<<" "<<f(Beta)<<std::endl;
  }


}//namespace
