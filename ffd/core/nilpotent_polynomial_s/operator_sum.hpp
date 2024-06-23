namespace ffd::nilpotent_polynomial_s {

template <class f, auto n, mode m>
auto& operator+=(P<f, n, m>& P0, P<f, n, m> const& P1) {
  using std::abs;
  if constexpr (m == unsafe) {
    for (ulong j = 0; j < (1ul << n); ++j) {
      P0[j] += P1[j];
    }  // for j in range(0, (1ul<<n))
  }
  if constexpr (m == safe) {
    auto const prec = std::max(P0.precision, P1.precision);
    P0.precision = prec;
    for (ulong j = 0; j < (1ul << n); ++j) {
      auto x = P0[j] + P1[j];
      auto y = abs(P0[j]) + abs(P1[j]);
      P0[j] = abs(x) < prec * y ? 0 : x;
    }  // for j in range(0, (1ul<<n))
  }

  return P0;
}

template <class f, auto n, mode m>
auto operator-(P<f, n, m> const& P0) {
  P<f, n, m> ret;
  ret.precision = P0.precision;
  for (ulong j = 0; j < (1ul << n); ++j) {
    ret[j] = -P0[j];
  }  // for j in range(0, (1ul<<n))
  return ret;
}

template <class f, auto n, mode m>
auto operator+(P<f, n, m> const& P0) {
  return P0;
}

template <class f, auto n, mode m>
auto& operator-=(P<f, n, m>& P0, P<f, n, m> const& P1) {
  return P0 += -P1;
}

template <class f, auto n, mode m>
auto operator+(P<f, n, m> const& P0, P<f, n, m> const& P1) {
  auto ret = P0;
  return ret += P1;
}

template <class f, auto n, mode m>
auto operator-(P<f, n, m> const& P0, P<f, n, m> const& P1) {
  auto ret = P0;
  return ret -= P1;
}

}  // namespace ffd::nilpotent_polynomial_s
