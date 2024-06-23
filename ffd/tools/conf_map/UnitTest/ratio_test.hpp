namespace ffd::conf_map::unit_test {

void ratio_test_test() {
  using namespace ffd::user_space;
  std::ifstream in_f = [] {
    std::stringstream ss;
    ss << "kx0_ky0_om0_i.series";
    return std::ifstream(ss.str());
  }();  // ifstream in_f
  ffd::taylor_polynomial::P<ENumber<Real>, 10> p;
  p = ENumber<Real>(0, 0);

  int order;
  Real val, err;
  while (in_f >> order >> val >> err) {
    p[order] = ENumber<Real>(val, err);
  }

  for (ulong j = 0; j <= p.order; ++j) {
    auto [val, err] = p[j].eval();
    std::cerr << j << " " << val << " " << err << "\n";
  }  // for j in range(0, p.order)

  auto const R0 = ratio_test(p, 6, 10);
  std::cerr << "R0 = " << R0 << "\n";
}

}  // namespace ffd::conf_map::unit_test
