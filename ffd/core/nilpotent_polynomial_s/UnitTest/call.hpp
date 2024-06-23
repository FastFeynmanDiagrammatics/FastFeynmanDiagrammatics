namespace ffd::nilpotent_polynomial_s::unit_test {

void UnitTest() {
  decl();
  test_multiply_divide<5, ffd::nilpotent_polynomial_s::safe>();
  test_multiply_divide<5, unsafe>();
}
}  // namespace ffd::nilpotent_polynomial_s::unit_test
