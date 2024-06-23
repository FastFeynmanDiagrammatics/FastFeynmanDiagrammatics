namespace ffd::alpha_f::unit_test {

template <uint ord, uint ord_>
void test_alpha() {
  ffd::array2d<Real, ord_, ord_> mat;
  ffd::array2d<ffd::truncated_polynomial::P<ord, Real>, ord_, ord_> mat_poly;
  for (uint i = 0; i < ord_; ++i) {
    for (uint j = 0; j < ord_; ++j) {
      mat(i, j) = (2. * ffd::user_space::Proba() - 1.) * (i != j || i >= ord);
      mat_poly(i, j) = mat(i, j);
      mat_poly(i, j)[1] = 1. * ((i == j) && i < ord);
    }
  }
  // ffd::math_tools::print_matrix(mat);

  auto minors = ffd::principal_minors::PrincipalMinorsBFS(mat);

  constexpr auto shift = (1 << ord_) - (1 << ord);
  auto alpha_minors = AlphaFunction<ord>(minors, shift);

  auto minors_poly = ffd::principal_minors::PrincipalMinorsBFS(mat_poly);

  for (uint j = 0; j < (1 << ord); ++j) {
    const uint card_j = __builtin_popcount(j);
    for (uint o = 0; o <= std::min(card_j, ord - card_j); ++o) {
      assert(std::abs(minors_poly[j + shift][o] - alpha_minors[j][o]) < 1e-10);
    }
    // std::cerr << j << " " << ToString(minors_poly[j + shift], 15) <<
    // std::endl
    //           << "  " << ToString(alpha_minors[j], 15) << std::endl;
  }
}

template <uint ord, uint ord_>
void test_signed_alpha() {
  std::array<Real, ord> sign;
  for (uint i = 0; i < ord; ++i) {
    sign[i] = ffd::random_distributions::RandomSign();
  }

  ffd::array2d<Real, ord_, ord_> mat;
  ffd::array2d<ffd::truncated_polynomial::P<ord, Real>, ord_, ord_> mat_poly;
  for (uint i = 0; i < ord_; ++i) {
    for (uint j = 0; j < ord_; ++j) {
      mat(i, j) = (2. * ffd::user_space::Proba() - 1.) * (i != j || i >= ord);
      mat_poly(i, j) = mat(i, j);
      mat_poly(i, j)[1] = sign[i] * ((i == j) && i < ord);
    }
  }
  // ffd::math_tools::print_matrix(mat);

  auto minors = ffd::principal_minors::PrincipalMinorsBFS(mat);
  constexpr auto shift = (1 << ord_) - (1 << ord);
  auto alpha_minors = SignedAlphaFunction<ord>(minors, sign, shift);

  auto minors_poly = ffd::principal_minors::PrincipalMinorsBFS(mat_poly);

  for (uint j = 0; j < (1 << ord); ++j) {
    const uint card_j = __builtin_popcount(j);
    for (uint o = 0; o <= std::min(card_j, ord - card_j); ++o) {
      assert(std::abs(minors_poly[j + shift][o] - alpha_minors[j][o]) < 1e-10);
    }
    // std::cerr << j << " " << ToString(minors_poly[j + shift], 10) <<
    // std::endl
    //           << "  " << ToString(alpha_minors[j], 10) << std::endl;
  }
}

template <uint ord, uint ord_>
void test_alpha_full() {
  ffd::array2d<Real, ord_, ord_> mat;
  ffd::array2d<ffd::truncated_polynomial::P<ord, Real>, ord_, ord_> mat_poly;
  for (uint i = 0; i < ord_; ++i) {
    for (uint j = 0; j < ord_; ++j) {
      mat(i, j) = (2. * ffd::user_space::Proba() - 1.) * (i != j || i >= ord);
      mat_poly(i, j) = mat(i, j);
      mat_poly(i, j)[1] = 1. * ((i == j) && i < ord);
    }
  }
  // ffd::math_tools::print_matrix(mat);

  auto minors = ffd::principal_minors::PrincipalMinorsBFS(mat);

  constexpr auto shift = (1 << ord_) - (1 << ord);
  auto alpha_minors = AlphaFunction<ord, true>(minors, shift);

  auto minors_poly = ffd::principal_minors::PrincipalMinorsBFS(mat_poly);

  for (uint j = 0; j < (1 << ord); ++j) {
    const uint card_j = __builtin_popcount(j);
    for (uint o = 0; o <= card_j; ++o) {
      assert(std::abs(minors_poly[j + shift][o] - alpha_minors[j][o]) < 1e-10);
    }
    // std::cerr << j << " " << ToString(minors_poly[j + shift], 8) << std::endl
    // << "  " << ToString(alpha_minors[j], 8) << std::endl;
  }
}

template <uint ord, uint ord_>
void test_signed_alpha_full() {
  std::array<Real, ord> sign;
  for (uint i = 0; i < ord; ++i) {
    sign[i] = ffd::random_distributions::RandomSign();
  }

  ffd::array2d<Real, ord_, ord_> mat;
  ffd::array2d<ffd::truncated_polynomial::P<ord, Real>, ord_, ord_> mat_poly;
  for (uint i = 0; i < ord_; ++i) {
    for (uint j = 0; j < ord_; ++j) {
      mat(i, j) = (2. * ffd::user_space::Proba() - 1.) * (i != j || i >= ord);
      mat_poly(i, j) = mat(i, j);
      mat_poly(i, j)[1] = sign[i] * ((i == j) && i < ord);
    }
  }
  // ffd::math_tools::print_matrix(mat);

  auto minors = ffd::principal_minors::PrincipalMinorsBFS(mat);
  constexpr auto shift = (1 << ord_) - (1 << ord);
  auto alpha_minors = SignedAlphaFunction<ord, true>(minors, sign, shift);

  auto minors_poly = ffd::principal_minors::PrincipalMinorsBFS(mat_poly);

  for (uint j = 0; j < (1 << ord); ++j) {
    const uint card_j = __builtin_popcount(j);
    for (uint o = 0; o <= card_j; ++o) {
      assert(std::abs(minors_poly[j + shift][o] - alpha_minors[j][o]) < 1e-10);
    }
    // std::cerr << j << " " << ToString(minors_poly[j + shift], 10) <<
    // std::endl
    //           << "  " << ToString(alpha_minors[j], 10) << std::endl;
  }
}

}  // namespace ffd::alpha_f::unit_test
