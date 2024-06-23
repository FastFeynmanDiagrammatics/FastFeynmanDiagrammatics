namespace ffd::taylor_polynomial {

template <class field, std::size_t n>
struct P {
  std::size_t static constexpr order = n;

  std::array<field, n + 1> coef;

  P() {}

  P(field x) {
    coef.fill(0);
    coef[0] = x;
  }

  explicit P(int j) {
    coef.fill(0);
    coef[0] = field(j);
  }

  inline field& operator[](std::size_t j) { return coef[j]; }

  inline field const& operator[](std::size_t j) const { return coef[j]; }

  field operator()(field x) const {
    field ret = 0.;
    field pow_x = 1.;
    for (ulong j = 0; j <= n; ++j) {
      ret += coef[j] * pow_x;
      pow_x *= x;
    }
    return ret;
  }

  auto& eval() { return coef; }
  auto& data() { return coef; }
};

template <class field, std::size_t n>
auto operator+(P<field, n> const& P1, P<field, n> const& P2) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] += P2[j];
  }
  return ret;
}

template <class field, std::size_t n>
auto& operator+=(P<field, n>& P1, P<field, n> const& P2) {
  P1 = P1 + P2;
  return P1;
}

template <class field, std::size_t n>
auto operator-(P<field, n> const& P1, P<field, n> const& P2) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] -= P2[j];
  }
  return ret;
}

template <class field, std::size_t n>
auto& operator-=(P<field, n>& P1, P<field, n> const& P2) {
  P1 = P1 - P2;
  return P1;
}

template <class field, std::size_t n>
auto operator*(P<field, n> const& P1, P<field, n> const& P2) {
  //  static_assert(n <= 10);
  std::size_t constexpr n_static_max = 10;
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
  if constexpr (n >= n_static_max) {
    ret.coef[10] = P1.coef[10] * P2.coef[0] + P1.coef[9] * P2.coef[1] +
                   P1.coef[8] * P2.coef[2] + P1.coef[7] * P2.coef[3] +
                   P1.coef[6] * P2.coef[4] + P1.coef[5] * P2.coef[5] +
                   P1.coef[4] * P2.coef[6] + P1.coef[3] * P2.coef[7] +
                   P1.coef[2] * P2.coef[8] + P1.coef[1] * P2.coef[9] +
                   P1.coef[0] * P2.coef[10];
  }
  if constexpr (n >= n_static_max + 1) {
    for (ulong j = n_static_max + 1; j <= n; ++j) {
      ret.coef[j] = 0;
      for (ulong k = 0; k < j + 1; ++k) {
        ret.coef[j] += P1.coef[k] * P2.coef[j - k];
      }  // for k in range(0, j+1)
    }    // for j in range(11, n+1)
  }

  return ret;
}

template <class field, std::size_t n>
auto& operator*=(P<field, n>& P1, P<field, n> const& P2) {
  P1 = P1 * P2;
  return P1;
}

template <class field, std::size_t n>
auto operator/(P<field, n> const& P1, field x) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] /= x;
  }
  return ret;
}

template <class field, std::size_t n>
auto& operator/=(P<field, n>& P1, field x) {
  P1 = P1 / x;
  return P1;
}

template <class field, std::size_t n>
auto operator/(field x, P<field, n> const& P1) {
  P<field, n> P2{x};
  auto ret = P2 / P1;
  return ret;
}

template <class field, std::size_t n>
auto operator/(P<field, n> const& P1, P<field, n> const& P3) {
  //  static_assert(n <= 10);
  std::size_t constexpr n_static_max = 10;
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
  if constexpr (n >= n_static_max) {
    ret.coef[10] -= P2.coef[10] * ret.coef[0] + P2.coef[9] * ret.coef[1] +
                    P2.coef[8] * ret.coef[2] + P2.coef[7] * ret.coef[3] +
                    P2.coef[6] * ret.coef[4] + P2.coef[5] * ret.coef[5] +
                    P2.coef[4] * ret.coef[6] + P2.coef[3] * ret.coef[7] +
                    P2.coef[2] * ret.coef[8] + P2.coef[1] * ret.coef[9];
  }
  if constexpr (n >= n_static_max + 1) {
    for (ulong j = n_static_max + 1; j <= n; ++j) {
      ret.coef[j] = 0;
      for (ulong k = 0; k < j; ++k) {
        ret.coef[j] -= ret.coef[k] * P2.coef[j - k];
      }  // for k in range(0, j)
    }    // for j in range(n_static_max+1, n)
  }

  return ret;
}

template <class field, std::size_t n>
auto& operator/=(P<field, n>& P1, P<field, n> const& P2) {
  P1 = P1 / P2;
  return P1;
}

template <class field, std::size_t n>
auto operator*(P<field, n> const& P1, field x) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] *= x;
  }
  return ret;
}

template <class field, std::size_t n>
auto operator*(field x, P<field, n> const& P1) {
  return P1 * x;
}

template <class field, std::size_t n>
auto& operator*=(P<field, n>& P1, field x) {
  P1 = P1 * x;
  return P1;
}

template <class field, std::size_t n>
auto operator+(P<field, n> const& P1, field x) {
  auto ret = P1;
  ret[0] += x;
  return ret;
}

template <class field, std::size_t n>
auto& operator+=(P<field, n>& P1, field x) {
  P1 = P1 + x;
  return P1;
}

template <class field, std::size_t n>
auto operator-(P<field, n> const& P1, field x) {
  auto ret = P1;
  ret[0] -= x;
  return ret;
}

template <class field, std::size_t n>
auto& operator-=(P<field, n>& P1, field x) {
  P1 = P1 - x;
  return P1;
}

template <class field, std::size_t n>
auto operator-(P<field, n> const& P1) {
  auto ret = P1;
  for (uint j = 0; j <= n; ++j) {
    ret[j] = -ret[j];
  }
  return ret;
}

template <class field, std::size_t n>
auto operator+(P<field, n> const& P1) {
  return P1;
}

template <class field, std::size_t n>
auto abs(P<field, n> const& P1) {
  using std::abs;
  return abs(P1[0]);
}

template <class field, std::size_t n>
auto pow(P<field, n> const& P1, int power) {
  if (power == 0) {
    return P<field, n>{field(1)};
  } else if (power < 0) {
    auto ret = field(1) / P1;
    for (ulong j = 0; j < std::abs(power) - 1; ++j) {
      ret /= P1;
    }  // for j in range(0, power)
    return ret;
  } else {
    auto ret = field(1) * P1;
    for (ulong j = 0; j < std::abs(power) - 1; ++j) {
      ret *= P1;
    }  // for j in range(0, power)
    return ret;
  }
}

}  // namespace ffd::taylor_polynomial
