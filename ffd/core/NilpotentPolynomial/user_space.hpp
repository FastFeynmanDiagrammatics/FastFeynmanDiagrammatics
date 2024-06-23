namespace ffd::user_space{

  using namespace ffd::SetT_wrapper;

  
  template<typename T, SetT SetType>
  using NilpotentPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<T, SetType>;

  using NilPolynomial = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;

  using NilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;

  using FullNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Full>;

  using EvenNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Even>;

  using OddNilPoly = ffd::nilpotent_polynomial::NilpotentPolynomial<Real, Odd>;



}//namespace ffd::user_space
