namespace ffd::nilpotent_polynomial_s {

template <class f, auto n, mode m>
auto operator/(P<f, n, m> const& P0, P<f, n, m> const& P1) {
  P<f, n, m> ret;
  f const one_over_Q0 = 1. / P1[0];
  for (ulong j = 0; j < (1ul << n); ++j) {
    ret[j] = P0[j] * one_over_Q0;
  }  // for j in range(0, (1ul<<n))

  auto const q = P1 * one_over_Q0;

  f const ret0 = ret[0];

  if constexpr (m == unsafe) {
    for (uint V = 1u; V < (1u << n); ++V) {
      ret[V] -= q[V] * ret0;
      uint const V_max = V - (1 << (31 - __builtin_clz(V)));
      for (uint S = ((V - 1) & V); S > V_max; S = ((S - 1) & V)) {
        uint const V_S = V - S;
        ret[V] -= q[S] * ret[V_S] + q[V_S] * ret[S];
      }
    }
  }

  if constexpr (m == safe) {
    auto const prec = std::max(P0.precision, P1.precision);
    ret.precision = prec;
    using std::abs;
    std::array<f, (1u << n)> ret_abs;
    for (uint V = 1u; V < (1u << n); ++V) {
      ret_abs[V] = abs(ret[V]);
      auto const x = q[V] * ret0;
      ret[V] -= x;
      ret_abs[V] += abs(x);
      uint const V_max = V - (1 << (31 - __builtin_clz(V)));
      for (uint S = ((V - 1) & V); S > V_max; S = ((S - 1) & V)) {
        uint const V_S = V - S;
        auto const y0 = q[S] * ret[V_S];
        auto const y1 = q[V_S] * ret[S];
        ret[V] -= y0 + y1;
        ret_abs[V] += abs(y0) + abs(y1);
      }
    }
    for (uint V = 1u; V < (1u << n); ++V) {
      if (abs(ret[V]) < prec * ret_abs[V]) {
        ret[V] = 0;
      }
    }
  }

  return ret;
}

template <class f, auto n, mode m>
auto& operator/=(P<f, n, m>& P0, P<f, n, m> const& P1) {
  return P0 = P0 / P1;
}

}  // namespace ffd::nilpotent_polynomial_s
