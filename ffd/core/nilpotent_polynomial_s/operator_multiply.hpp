namespace ffd::nilpotent_polynomial_s {

template <class f, auto n, mode m>
auto operator*(P<f, n, m> const& P0, P<f, n, m> const& P1) {
  P<f, n, m> ret;

  if constexpr (m == unsafe) {
    ret[0] = P0[0] * P1[0];
    for (uint V = 1u; V < (1ul << n); ++V) {
      ret[V] = 0;
      uint const V_max = V - (1 << (31 - __builtin_clz(V)));
      for (uint S = V; S > V_max; S = ((S - 1) & V)) {
        uint const V_S = V - S;
        ret[V] += P0[S] * P1[V_S] + P0[V_S] * P1[S];
      }
    }
  }

  if constexpr (m == safe) {
    auto const prec = std::max(P0.precision, P1.precision);
    ret.precision = prec;
    using std::abs;
    ret[0] = P0[0] * P1[0];
    for (uint V = 1u; V < (1ul << n); ++V) {
      ret[V] = 0;
      auto ret_abs = abs(ret[V]);
      uint const V_max = V - (1 << (31 - __builtin_clz(V)));
      for (uint S = V; S > V_max; S = ((S - 1) & V)) {
        uint const V_S = V - S;
        auto x0 = P0[S] * P1[V_S];
        auto x1 = P0[V_S] * P1[S];
        ret[V] += x0 + x1;
        ret_abs += abs(x0) + abs(x1);
      }
      if (abs(ret[V]) < prec * ret_abs)
        ret[V] = 0;
    }
  }

  return ret;
}

template <class f, auto n, mode m>
auto& operator*=(P<f, n, m>& P0, P<f, n, m> const& P1) {
  return P0 = P0 * P1;
}

}  // namespace ffd::nilpotent_polynomial_s
