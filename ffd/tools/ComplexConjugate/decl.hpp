

namespace ffd::complex_conjugate{

  template<typename T>
  inline T ComplexConjugate(T x_){return x_;}

  template<typename T>
  inline std::complex<T> ComplexConjugate(std::complex<T> x_){return std::conj(x_);}
  
}//namespace
