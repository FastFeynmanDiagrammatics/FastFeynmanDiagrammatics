namespace ffd::truncated_polynomial {

template <uint n, class field_d = Real>
class P {
 public:
  std::array<field_d, n + 1> coef;

  P() { coef.fill(0.); }  // maybe change

  P(field_d x) {
    coef.fill(0);
    coef[0] = x;
  }

  inline field_d& operator[](std::size_t j) { return coef[j]; }

  inline field_d const& operator[](std::size_t j) const { return coef[j]; }

  field_d operator()(field_d x) const {  // should be done better
    field_d ret = 0.;
    field_d pow_x = 1.;
    for (uint j = 0; j <= n; ++j) {
      ret += coef[j] * pow_x;
      pow_x *= x;
    }
  }
};

template <uint n, class field_d>
auto operator+(P<n, field_d> const& P1, P<n, field_d> const& P2) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] += P2[j];
  }
  return ret;
}

template <uint n, class field_d>
auto& operator+=(P<n, field_d>& P1, P<n, field_d> const& P2) {
  P1 = P1 + P2;
  return P1;
}

template <uint n, class field_d>
auto operator-(P<n, field_d> const& P1, P<n, field_d> const& P2) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] -= P2[j];
  }
  return ret;
}

template <uint n, class field_d>
auto& operator-=(P<n, field_d>& P1, P<n, field_d> const& P2) {
  P1 = P1 - P2;
  return P1;
}

template <uint n, class field_d>
auto operator*(P<n, field_d> const& P1, P<n, field_d> const& P2) {
  static_assert(n <= 10);
  auto ret = P1;

  if constexpr (n >= 0) {
    ret.coef[0] = P1.coef[0] * P2.coef[0];
  }
  if constexpr (n >= 1) {
    ret.coef[1] = P1.coef[1] * P2.coef[0] + P1.coef[0] * P2.coef[1];
  }
  if constexpr (n >= 2) {
    ret.coef[2] = P1.coef[2] * P2.coef[0] + P1.coef[1] * P2.coef[1] +
                  P1.coef[0] * P2.coef[2];
  }
  if constexpr (n >= 3) {
    ret.coef[3] = P1.coef[3] * P2.coef[0] + P1.coef[2] * P2.coef[1] +
                  P1.coef[1] * P2.coef[2] + P1.coef[0] * P2.coef[3];
  }
  if constexpr (n >= 4) {
    ret.coef[4] = P1.coef[4] * P2.coef[0] + P1.coef[3] * P2.coef[1] +
                  P1.coef[2] * P2.coef[2] + P1.coef[1] * P2.coef[3] +
                  P1.coef[0] * P2.coef[4];
  }
  if constexpr (n >= 5) {
    ret.coef[5] = P1.coef[5] * P2.coef[0] + P1.coef[4] * P2.coef[1] +
                  P1.coef[3] * P2.coef[2] + P1.coef[2] * P2.coef[3] +
                  P1.coef[1] * P2.coef[4] + P1.coef[0] * P2.coef[5];
  }
  if constexpr (n >= 6) {
    ret.coef[6] = P1.coef[6] * P2.coef[0] + P1.coef[5] * P2.coef[1] +
                  P1.coef[4] * P2.coef[2] + P1.coef[3] * P2.coef[3] +
                  P1.coef[2] * P2.coef[4] + P1.coef[1] * P2.coef[5] +
                  P1.coef[0] * P2.coef[6];
  }
  if constexpr (n >= 7) {
    ret.coef[7] = P1.coef[7] * P2.coef[0] + P1.coef[6] * P2.coef[1] +
                  P1.coef[5] * P2.coef[2] + P1.coef[4] * P2.coef[3] +
                  P1.coef[3] * P2.coef[4] + P1.coef[2] * P2.coef[5] +
                  P1.coef[1] * P2.coef[6] + P1.coef[0] * P2.coef[7];
  }
  if constexpr (n >= 8) {
    ret.coef[8] = P1.coef[8] * P2.coef[0] + P1.coef[7] * P2.coef[1] +
                  P1.coef[6] * P2.coef[2] + P1.coef[5] * P2.coef[3] +
                  P1.coef[4] * P2.coef[4] + P1.coef[3] * P2.coef[5] +
                  P1.coef[2] * P2.coef[6] + P1.coef[1] * P2.coef[7] +
                  P1.coef[0] * P2.coef[8];
  }
  if constexpr (n >= 9) {
    ret.coef[9] = P1.coef[9] * P2.coef[0] + P1.coef[8] * P2.coef[1] +
                  P1.coef[7] * P2.coef[2] + P1.coef[6] * P2.coef[3] +
                  P1.coef[5] * P2.coef[4] + P1.coef[4] * P2.coef[5] +
                  P1.coef[3] * P2.coef[6] + P1.coef[2] * P2.coef[7] +
                  P1.coef[1] * P2.coef[8] + P1.coef[0] * P2.coef[9];
  }
  if constexpr (n >= 10) {
    ret.coef[10] = P1.coef[10] * P2.coef[0] + P1.coef[9] * P2.coef[1] +
                   P1.coef[8] * P2.coef[2] + P1.coef[7] * P2.coef[3] +
                   P1.coef[6] * P2.coef[4] + P1.coef[5] * P2.coef[5] +
                   P1.coef[4] * P2.coef[6] + P1.coef[3] * P2.coef[7] +
                   P1.coef[2] * P2.coef[8] + P1.coef[1] * P2.coef[9] +
                   P1.coef[0] * P2.coef[10];
  }

  return ret;
}

template <uint n, class field_d>
auto& operator*=(P<n, field_d>& P1, P<n, field_d> const& P2) {
  P1 = P1 * P2;
  return P1;
}

template <uint n, class field_d>
auto operator/(P<n, field_d> const& P1, field_d x) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] /= x;
  }
  return ret;
}

template <uint n, class field_d>
auto& operator/=(P<n, field_d>& P1, field_d x) {
  P1 = P1 / x;
  return P1;
}

template <uint n, class field_d>
auto operator/(field_d x, P<n, field_d> const& P1) {
  P<n, field_d> P2 = x;
  auto ret = P2 / P1;
  return ret;
}

template <uint n, class field_d>
auto operator/(P<n, field_d> const& P1, P<n, field_d> const& P3) {
  static_assert(n <= 10);
  auto ret = P1 / P3.coef[0];
  auto const P2 = P3 / P3.coef[0];

  if constexpr (n >= 1) {
    ret.coef[1] -= P2.coef[1] * ret.coef[0];
  }
  if constexpr (n >= 2) {
    ret.coef[2] -= P2.coef[2] * ret.coef[0] + P2.coef[1] * ret.coef[1];
  }
  if constexpr (n >= 3) {
    ret.coef[3] -= P2.coef[3] * ret.coef[0] + P2.coef[2] * ret.coef[1] +
                   P2.coef[1] * ret.coef[2];
  }
  if constexpr (n >= 4) {
    ret.coef[4] -= P2.coef[4] * ret.coef[0] + P2.coef[3] * ret.coef[1] +
                   P2.coef[2] * ret.coef[2] + P2.coef[1] * ret.coef[3];
  }
  if constexpr (n >= 5) {
    ret.coef[5] -= P2.coef[5] * ret.coef[0] + P2.coef[4] * ret.coef[1] +
                   P2.coef[3] * ret.coef[2] + P2.coef[2] * ret.coef[3] +
                   P2.coef[1] * ret.coef[4];
  }
  if constexpr (n >= 6) {
    ret.coef[6] -= P2.coef[6] * ret.coef[0] + P2.coef[5] * ret.coef[1] +
                   P2.coef[4] * ret.coef[2] + P2.coef[3] * ret.coef[3] +
                   P2.coef[2] * ret.coef[4] + P2.coef[1] * ret.coef[5];
  }
  if constexpr (n >= 7) {
    ret.coef[7] -= P2.coef[7] * ret.coef[0] + P2.coef[6] * ret.coef[1] +
                   P2.coef[5] * ret.coef[2] + P2.coef[4] * ret.coef[3] +
                   P2.coef[3] * ret.coef[4] + P2.coef[2] * ret.coef[5] +
                   P2.coef[1] * ret.coef[6];
  }
  if constexpr (n >= 8) {
    ret.coef[8] -= P2.coef[8] * ret.coef[0] + P2.coef[7] * ret.coef[1] +
                   P2.coef[6] * ret.coef[2] + P2.coef[5] * ret.coef[3] +
                   P2.coef[4] * ret.coef[4] + P2.coef[3] * ret.coef[5] +
                   P2.coef[2] * ret.coef[6] + P2.coef[1] * ret.coef[7];
  }
  if constexpr (n >= 9) {
    ret.coef[9] -= P2.coef[9] * ret.coef[0] + P2.coef[8] * ret.coef[1] +
                   P2.coef[7] * ret.coef[2] + P2.coef[6] * ret.coef[3] +
                   P2.coef[5] * ret.coef[4] + P2.coef[4] * ret.coef[5] +
                   P2.coef[3] * ret.coef[6] + P2.coef[2] * ret.coef[7] +
                   P2.coef[1] * ret.coef[8];
  }
  if constexpr (n >= 10) {
    ret.coef[10] -= P2.coef[10] * ret.coef[0] + P2.coef[9] * ret.coef[1] +
                    P2.coef[8] * ret.coef[2] + P2.coef[7] * ret.coef[3] +
                    P2.coef[6] * ret.coef[4] + P2.coef[5] * ret.coef[5] +
                    P2.coef[4] * ret.coef[6] + P2.coef[3] * ret.coef[7] +
                    P2.coef[2] * ret.coef[8] + P2.coef[1] * ret.coef[9];
  }

  return ret;
}

template <uint n, class field_d>
auto& operator/=(P<n, field_d>& P1, P<n, field_d> const& P2) {
  P1 = P1 / P2;
  return P1;
}

template <uint n, class field_d>
auto operator*(P<n, field_d> const& P1, field_d x) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] *= x;
  }
  return ret;
}

template <uint n, class field_d>
auto operator*(field_d x, P<n, field_d> const& P1) {
  return P1 * x;
}

template <uint n, class field_d>
auto& operator*=(P<n, field_d>& P1, field_d x) {
  P1 = P1 * x;
  return P1;
}

template <uint n, class field_d>
auto operator+(P<n, field_d> const& P1, field_d x) {
  auto ret = P1;
  ret[0] += x;
  return ret;
}

template <uint n, class field_d>
auto& operator+=(P<n, field_d>& P1, field_d x) {
  P1 = P1 + x;
  return P1;
}

template <uint n, class field_d>
auto operator-(P<n, field_d> const& P1, field_d x) {
  auto ret = P1;
  ret[0] -= x;
  return ret;
}

template <uint n, class field_d>
auto& operator-=(P<n, field_d>& P1, field_d x) {
  P1 = P1 - x;
  return P1;
}

template <uint n, class field_d>
auto operator-(P<n, field_d> const& P1) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] = -ret[j];
  }
  return ret;
}

template <uint n, class field_d>
auto operator+(P<n, field_d> const& P1) {
  return P1;
}

template <uint n, class field_d>
auto operator<(P<n, field_d> const& P1, P<n, field_d> const& P2) {
  return (P1[0] < P2[0]);
}

template <uint n, class field_d>
auto operator>(P<n, field_d> const& P1, P<n, field_d> const& P2) {
  return (P1[0] > P2[0]);
}

template <uint n, class field_d>
auto operator>(P<n, field_d> const& P1, field_d x) {
  return (P1[0] > x);
}

template <uint n, class field_d>
auto operator<(P<n, field_d> const& P1, field_d x) {
  return (P1[0] < x);
}

template <uint n, class field_d>
auto abs(P<n, field_d> const& P1) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] = std::abs(ret[j]);
  }
  return ret;
}

template <uint n, class field_d>
auto shift(P<n, field_d> const& P0, field_d x0) {
  using std::exp, std::pow, std::lgamma;
  P<n, field_d> ret;
  ret.coef.fill(0);
  for (int j = 0; j < n + 1; ++j) {
    field_d pow_x0 = 1;
    for (int k = 0; k < j + 1; ++k, pow_x0 *= x0) {
      ret.coef[j - k] +=
          exp(lgamma(j + 1.) - lgamma(k + 1.) - lgamma(j - k + 1.)) * pow_x0 *
          P0.coef[j];
    }  // for k in range(0, j+1)
  }    // for j in range(0, n+1)
  return ret;
}

}  // namespace ffd::truncated_polynomial
