namespace ffd::user_space{

  template<typename T>
  using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<T, Full>;

  using NilPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;

  using NilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;

  using FullNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;

  using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;

  using OddNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Odd>;

}//namespace ffd::user_space
