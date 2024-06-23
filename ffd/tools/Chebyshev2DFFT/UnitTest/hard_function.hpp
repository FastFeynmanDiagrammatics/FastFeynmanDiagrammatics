

namespace ffd::chebyshev_2d_fft::unit_test{

  struct hard_function_x{
    
    Real operator()(Real x, Real y){
      return sin(y) + x*sin(69*sin(72*x*x) + 12);
    }
    
  };


}//namespace
