namespace ffd::inverse_matrix_gauss {

template <typename T>
bool NotZero(T const& x) {
  using std::abs;
  return abs(x) > 100 * std::numeric_limits<Real>::min();
}

bool NotZero(ffd::nilpotent_polynomial::NilpotentPolynomial<Real> const& x) {
  return x.NotZero;
}

}  // namespace ffd::inverse_matrix_gauss
