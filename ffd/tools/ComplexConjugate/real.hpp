

namespace ffd::real_complex_tools{

  template<typename T>
  inline T real(T x_){return x_;}

  template<typename T>
  inline T real(std::complex<T> x_){return std::real(x_);}



}//namespace
