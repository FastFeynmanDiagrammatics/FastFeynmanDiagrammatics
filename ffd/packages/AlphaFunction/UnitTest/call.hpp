namespace ffd::alpha_f::unit_test {

void UnitTest() {
  test_alpha<3, 3>();
  test_alpha<3, 4>();
  test_alpha<3, 5>();
  test_alpha<6, 7>();

  test_alpha_full<3, 3>();
  test_alpha_full<3, 4>();
  test_alpha_full<3, 5>();
  // test_alpha_full<6, 7>();

  test_signed_alpha<3, 3>();
  test_signed_alpha<3, 4>();
  test_signed_alpha<3, 5>();
  test_signed_alpha<6, 7>();

  test_signed_alpha_full<3, 3>();
  test_signed_alpha_full<3, 4>();
  test_signed_alpha_full<3, 5>();
  // test_signed_alpha_full<6, 7>();
}

}  // namespace ffd::alpha_f::unit_test
