namespace ffd::conf_map {

template <class poly_t>
auto ratio_test(poly_t const& P, int j_min, int j_max) {
  using std::abs, std::exp, std::log;
  using field = std::decay_t<decltype(P[0].eval().first)>;
  field Y = 0, X = 0, Z = 0, I = 0, W = 0;
  for (ulong j = j_min; j < j_max + 1; ++j) {
    auto [val, err] = P[j].eval();
    auto Err = err / abs(val);
    auto w = 1 / (Err * Err);
    Y += log(abs(val)) * w;
    X += j * w;
    Z += j * log(abs(val)) * w;
    I += j * j * w;
    W += w;
  }  // for j in range(j_min, j_max+1)
  //  std::cerr << Z / W << " " << X / W << " " << Y / W << " " << I / W
  //        << " aaa\n";
  return exp(-(Z / W - X * Y / W / W) / (I / W - X * X / W / W));
}

}  // namespace ffd::conf_map
