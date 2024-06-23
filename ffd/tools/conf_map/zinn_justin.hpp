namespace ffd::conf_map {

auto zinn_justin_g(Real R) {
  auto ret = [R]<class field>(field w) {
    using std::abs, std::pow;
    return field(4 * R) * w * pow(field(1) - w, -2);
  };
  return ret;
}

}  // namespace ffd::conf_map
