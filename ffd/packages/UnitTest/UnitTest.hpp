#ifdef FFD_UNIT_TEST_FLAG

namespace ffd::packages_unit_test {

int UnitTest() {
  ffd::user_space::Timer packages_unit_test_time;
  packages_unit_test_time.ini();

  ffd::counter_unit_test = 0;

  std::cerr << "test# = ";

  unit_test_runner(ffd::packages_math::unit_test::UnitTest);

  unit_test_runner(ffd::relation_sets::unit_test::UnitTest);

  unit_test_runner(ffd::class_tuple::unit_test::UnitTest);

  //    unit_test_runner(ffd::periodic_coordinate::unit_test::UnitTest);

  unit_test_runner(ffd::lattice::unit_test::UnitTest);

  unit_test_runner(ffd::imaginary_time::unit_test::UnitTest);

  unit_test_runner(ffd::imaginary_time_lattice::unit_test::UnitTest);

  unit_test_runner(ffd::user_space::quadratic_action::unit_test::UnitTest);

  unit_test_runner(ffd::wick_block_oracle::unit_test::UnitTest);

  unit_test_runner(ffd::imaginary_time_convolution::unit_test::UnitTest);

  unit_test_runner(ffd::imaginary_time_convolution_2d::unit_test::UnitTest);

  //    unit_test_runner(ffd::rpa_ladder::unit_test::UnitTest);

  unit_test_runner(
      ffd::user_space::imaginary_time_lattice_propagator::unit_test::UnitTest);

  unit_test_runner(ffd::imaginary_time_lattice_proposer::unit_test::UnitTest);

  std::cerr << " ||  ";

  return int(packages_unit_test_time.elapsed() + .5);
}

}  // namespace ffd::packages_unit_test

#endif
