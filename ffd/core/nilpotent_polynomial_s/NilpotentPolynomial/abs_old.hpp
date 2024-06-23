namespace ffd::nilpotent_polynomial{

  template<typename T>
  auto abs(NilpotentPolynomial<T> const& P) noexcept{
    using std::abs;
    return abs(P[0]);
  }

}//namespace ffd::nilpotent_polynomial
