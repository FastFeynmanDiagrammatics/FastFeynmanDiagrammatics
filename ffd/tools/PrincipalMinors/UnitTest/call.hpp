namespace ffd::principal_minors::unit_test {

void UnitTest() {
  matrix_4_by_4_Leading();

  matrix_6_by_6_Real();

  test_bounded();

  matrix_3_by_3_Real();

  matrix_3_by_3_classic();

  matrix_4_by_4_Real();

  matrix_6_by_6_Real();

  matrix_3_by_3_Poly();

  matrix_4_by_4_step();

  matrix_4_by_4_step_0_diag();

  matrix_6_by_6_step();

  matrix_6_by_6_step_2();

  matrix_16_by_16_step_2();

  test_max_order();
}

}  // namespace ffd::principal_minors::unit_test
