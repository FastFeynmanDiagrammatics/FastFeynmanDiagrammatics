

namespace ffd::chebyshev_fft::unit_test{

  Real func_exp_sin_test(Real x){
    return exp(-3*x)*sin(10*x)/(1+x*x);
  }

  
  void
  instantiation(){
    ChebyshevFFT P;

    ChebyshevFFT Q = 1.;

    ChebyshevFFT R(func_exp_sin_test, {-.3, 1.});

    // std::cerr<<std::setprecision(16)<<R(.3)<<" "<<func_exp_sin_test(.3)<<std::endl;
  }


}//namespace
