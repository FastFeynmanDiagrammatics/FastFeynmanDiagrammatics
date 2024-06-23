namespace ffd::imaginary_time_convolution::unit_test{
  
  void test_convolution_two_G0(Real beta, Real E1, Real E2, bool is_fermion, Real precision){
    using std::abs;

    
    auto f_E = [beta, is_fermion](Real tau, Real E)->Real{
		 return -exp(-tau*E)/(1+(-1+2*is_fermion)*exp(-beta*E));};
    auto f1 = [f_E, E1](Real tau)->Real{return f_E(tau, E1);};
    auto f2 = [f_E, E2](Real tau)->Real{return f_E(tau, E2);};
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Real> P1(f1, {0, beta}, precision);
    ffd::chebyshev_polynomial::ChebyshevPolynomial<Real> P2(f2, {0, beta}, precision);

    //((w-E1)(w-E2))^{-1}=(E1-E2)^{-1}((w-E1)^{-1}-(w-E2)^{-1})
    auto P3 = Convolute1And2FromZeroTo3<Real>(P1, P2, beta, is_fermion, precision);
    auto f3 = [f1, f2, E1, E2](Real tau)->Real{return (1./(E1-E2))*(f1(tau)-f2(tau));};

    const int N_test = 20;
    for(int j=0; j<N_test; ++j){
      Real tau = (j+.5)*beta/N_test;
      assert( abs( f3(tau) - P3(tau) ) < precision );
    }
    
  }
  
}//namespace
