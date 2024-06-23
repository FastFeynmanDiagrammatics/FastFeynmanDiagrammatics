

namespace ffd::imaginary_time_convolution::unit_test{

  void UnitTest(){
    using namespace std;
    
    Real a = 23.12312, b = 2.123382349;
    Real beta = 2.;
    bool is_fermion;
    Real precision = 1e-12;
    auto f = [a](Real){return a;};
    auto g = [b](Real){return b;};

    auto f_g = Convolute1And2FromZeroTo3<Real>(f, g, beta, is_fermion=false, precision);

    Real tau = .41234123;

    assert( abs( f_g(tau) - a*b*beta) < precision );


    
    Real Energy1=1.2312, Energy2=2.321123;
    test_convolution_two_G0(Energy1, Energy2, beta, is_fermion = true, precision);

    

    
    test_convolution_two_G0(Energy1 = 0.123124, Energy2 = 0.432312312, beta, is_fermion = false, precision);
    
    

  }
    
}//namespace
