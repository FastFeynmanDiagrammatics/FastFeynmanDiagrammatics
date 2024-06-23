namespace ffd::conf_map::unit_test {

void map_poly_test() {
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

  auto zj = zinn_justin_g(R0);

  auto pw = map_poly(zj, p);

  auto U = 10.;
  auto w = root(zj, U);
  std::cerr << "w(" << U << ")=" << w << "\n";

  for (ulong j = 0; j <= pw.order; ++j) {
    auto [val, err] = pw[j].eval();
    std::cerr << j << " " << val * std::pow(w, j) << " " << err * std::pow(w, j)
              << "\n";
  }  // for j in range(0, p.order)

  std::vector<Real> PW(11);
  for (ulong j = 0; j < pw.order + 1; ++j) {
    PW[j] = pw[j].value;
  }  // for j in range(0, pw.order+1)
  auto [pp, qq] = ffd::pade::eval(4, 4, PW);
  auto const dU = 0.25;
  for (ulong u = 0; u < 48; ++u) {
    auto const U_ = u * dU;
    auto w_U = root(zj, U_);
    auto const res = ffd::pade::eval_rational_function(pp, qq, w_U);
    std::cerr << "w=" << w_U << ", U=" << U_ << ", " << res << "\n";
  }
}

}  // namespace ffd::conf_map::unit_test
