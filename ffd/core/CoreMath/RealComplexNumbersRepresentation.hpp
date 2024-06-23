namespace ffd{

  using Real = double;

  using Complex = std::complex<Real>;

}


namespace ffd::user_space{

  using ffd::Real;

  using ffd::Complex;


  constexpr Complex I = Complex(0., 1.);

  
  Complex expI(Real x){
    return std::exp(Complex(0., x));
  }
  

}
