namespace ffd::determinant::unit_test {

void UnitTest() {
  matrix2();

  matrix3();

  matrix4();

  matrix5();

  matrix6();

  matrix7();

  // UnitTest/NilpotentPolynomial.hpp
  two_by_two();

  three_by_three();

  four_by_four();

  // five_by_five(); //too many floating point operations,
  // exact result often unrealiable

  compare_gauss_pp_cp<30>();
  compare_gauss_quad<50>();
}

}  // namespace ffd::determinant::unit_test
