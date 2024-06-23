namespace ffd::linear_cubic_interpolation {

template <typename element_t, typename array_t>
auto Interpolate(const array_t& Axyz,
                 const element_t& val,
                 const element_t& step) {
  const int px = val[0] / step[0];
  const int py = val[1] / step[1];
  const int pz = val[2] / step[2];
  const auto px1 = px + 1;
  const auto py1 = py + 1;
  const auto pz1 = pz + 1;
  const auto x = val[0] - px;
  const auto y = val[1] - py;
  const auto z = val[2] - pz;
  const auto v1 = Axyz(px, py, pz);
  const auto v2 = Axyz(px1, py, pz);
  const auto v3 = Axyz(px, py1, pz);
  const auto v4 = Axyz(px1, py1, pz);
  const auto v5 = Axyz(px, py, pz1);
  const auto v6 = Axyz(px1, py, pz1);
  const auto v7 = Axyz(px, py1, pz1);
  const auto v8 = Axyz(px1, py1, pz1);
  const auto x12 = v1 * (1 - x) + v2 * x;
  const auto x34 = v3 * (1 - x) + v4 * x;
  const auto x56 = v5 * (1 - x) + v6 * x;
  const auto x78 = v7 * (1 - x) + v8 * x;
  const auto y1234 = x12 * (1 - y) + x34 * y;
  const auto y5678 = x56 * (1 - y) + x78 * y;
  return (y1234 * (1 - z) + y5678 * z);
}

}  // namespace ffd::linear_cubic_interpolation
