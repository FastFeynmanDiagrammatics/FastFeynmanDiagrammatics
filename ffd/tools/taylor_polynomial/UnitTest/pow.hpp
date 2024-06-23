namespace ffd::taylor_polynomial::unit_test {

void pow() {
  P<Real, 10> p;
  for (ulong j = 0; j <= p.order; ++j) {
    p[j] = 1;
  }  // for j in range(0, )

  auto p0 = pow(p, 0);
  auto p1 = pow(p, 1);
  auto p2 = pow(p, 2);
  auto p3 = pow(p, 3);
  auto p_1 = pow(p, -1);
  auto p_2 = pow(p, -2);

  for (ulong j = 0; j < p.order; ++j) {
    std::cerr << j << " p=" << p[j] << " p0=" << p0[j] << " p1=" << p1[j]
              << " p2=" << p2[j] << " p3=" << p3[j] << " p_1=" << p_1[j]
              << " p_2=" << p_2[j] << "\n";
  }  // for j in range(0, p.order)
}

}  // namespace ffd::taylor_polynomial::unit_test
