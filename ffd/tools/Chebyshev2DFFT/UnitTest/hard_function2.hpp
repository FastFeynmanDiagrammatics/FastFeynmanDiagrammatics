

namespace ffd::chebyshev_2d_fft::unit_test{

  struct hard_function2{

    Real
    operator()(Real x, Real y) const{
      std::swap(x, y);
      return tanh(x*sqrt(y)) + sin(x-y*y/(1+x*x))*exp(-exp(x*y));
    }

  };


}//namespace
