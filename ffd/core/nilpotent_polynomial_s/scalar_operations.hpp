namespace ffd::nilpotent_polynomial_s {

template <class f, auto n, mode m>
auto& operator+=(P<f, n, m>& P0, f x) {
  using std::abs;
  if constexpr (m == unsafe) {
    P0[0] += x;
  }
  if constexpr (m == safe) {
    auto y = x + P0[0];
    auto z = abs(x) + abs(P0[0]);
    P0[0] = abs(y) < P0.precision * z ? 0 : y;
  }
  return P0;
}  // namespace ffd::nilpotent_polynomial_s

template <class f, auto n, mode m>
auto& operator-=(P<f, n, m>& P0, f x) {
  return P0 += -x;
}

template <class f, auto n, mode m>
auto operator+(P<f, n, m> const& P0, f x) {
  auto ret = P0;
  return ret += x;
}

template <class f, auto n, mode m>
auto operator-(P<f, n, m> const& P0, f x) {
  auto ret = P0;
  return ret -= x;
}

template <class f, auto n, mode m>
auto& operator*=(P<f, n, m>& P0, f x) {
  using std::abs;
  for (ulong j = 0; j < (1ul << n); ++j) {
    P0[j] *= x;
  }  // for j in range(0, (1ul<<n))
  return P0;
}

template <class f, auto n, mode m>
auto& operator/=(P<f, n, m>& P0, f x) {
  using std::abs;
  for (ulong j = 0; j < (1ul << n); ++j) {
    P0[j] /= x;
  }  // for j in range(0, (1ul<<n))
  return P0;
}

template <class f, auto n, mode m>
auto operator*(P<f, n, m> const& P0, f x) {
  auto ret = P0;
  return ret *= x;
}

template <class f, auto n, mode m>
auto operator/(P<f, n, m> const& P0, f x) {
  auto ret = P0;
  return ret /= x;
}

}  // namespace ffd::nilpotent_polynomial_s
