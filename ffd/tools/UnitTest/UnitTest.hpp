namespace ffd::tools_unit_test {

int UnitTest() {
  ffd::user_space::Timer tools_unit_test_time;
  tools_unit_test_time.ini();

  ffd::counter_unit_test = 0;

  std::cerr << "test# = ";

  unit_test_runner(ffd::math_tools::unit_test::UnitTest);

  unit_test_runner(ffd::chebyshev_polynomial::unit_test::UnitTest);

  unit_test_runner(ffd::integrate_clenshaw_curtis::unit_test::UnitTest);

  unit_test_runner(ffd::integral_function::unit_test::UnitTest);

  unit_test_runner(ffd::chebyshev_polynomial::unit_test::UnitTest);

  unit_test_runner(ffd::eigenvalues_vectors::unit_test::UnitTest);

  unit_test_runner(ffd::inverse_matrix::unit_test::UnitTest);

  unit_test_runner(ffd::vector_range::unit_test::UnitTest);

  unit_test_runner(ffd::cartesian_product_array::unit_test::UnitTest);

  unit_test_runner(ffd::combination::unit_test::UnitTest);

  unit_test_runner(ffd::cartesian_product::unit_test::UnitTest);

  unit_test_runner(ffd::fft::unit_test::UnitTest);

  unit_test_runner(ffd::chebyshev_fft::unit_test::UnitTest);

  unit_test_runner(ffd::chebyshev_2d_fft::unit_test::UnitTest);

  unit_test_runner(ffd::find_root::unit_test::UnitTest);

  unit_test_runner(ffd::integrate::unit_test::UnitTest);

  unit_test_runner(ffd::contour_derivative::unit_test::UnitTest);

  unit_test_runner(ffd::diagmc_sampling::unit_test::UnitTest);

  unit_test_runner(ffd::enumerate::unit_test::UnitTest);

  unit_test_runner(ffd::user_space::safe_array::unit_test::UnitTest);

  unit_test_runner(ffd::user_space::safe_map::unit_test::UnitTest);

  unit_test_runner(ffd::gauss_pfaffian::unit_test::UnitTest);

  unit_test_runner(ffd::principal_minors::unit_test::UnitTest);

  unit_test_runner(ffd::determinant::unit_test::UnitTest);

  unit_test_runner(ffd::gauss_kronrod::unit_test::UnitTest);

  unit_test_runner(ffd::inverse_matrix_gauss::unit_test::UnitTest);

  std::cerr << " ||  ";

  return int(tools_unit_test_time.elapsed() + .5);
}

}  // namespace ffd::tools_unit_test
