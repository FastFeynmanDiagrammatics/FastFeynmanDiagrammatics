namespace ffd::truncated_polynomial {

template <uint n, class field_d>
auto ToString(P<n, field_d> const& P, int precision = 10) {
  std::stringstream ss;
  // for(uint j=0; j<=P.order; ++j){
  for (uint j = 0; j <= n; ++j) {
    ss << ' ' << std::setprecision(precision) << std::showpos << P[j] << "z^"
       << j;
  }
  return ss.str();
}

}  // namespace ffd::truncated_polynomial
