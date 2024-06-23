namespace ffd::round_off {

template <class field>
struct R {
  field first;
  field second;
  R() = default;
  R(field x) : first(x) {
    second = std::numeric_limits<field>::epsilon() * (x > 0 ? x : -x);
  }
  R(field x, field y) : first(x), second(y) {}
};

template <class f>
inline auto operator+(R<f> x, R<f> y) {
  return R<f>{x.first + y.first, x.second + y.second};
}

template <class f>
inline auto operator+(R<f> x, f y) {
  return R<f>{x.first + y, x.second};
}

template <class f>
inline auto operator+(f x, R<f> y) {
  return R<f>{x + y.first, y.second};
}

template <class f>
inline auto& operator+=(R<f>& x, R<f> y) {
  x.first += y.first;
  x.second += y.second;
  return x;
}

template <class f>
inline auto& operator+=(R<f>& x, f y) {
  x.first += y;
  return x;
}

template <class f>
inline auto operator-(R<f> x, R<f> y) {
  return R<f>{x.first - y.first, x.second + y.second};
}

template <class f>
inline auto operator-(R<f> x, f y) {
  return R<f>{x.first - y, x.second};
}

template <class f>
inline auto operator-(f x, R<f> y) {
  return R<f>{x - y.first, y.second};
}

template <class f>
inline auto& operator-=(R<f>& x, R<f> y) {
  x.first -= y.first;
  x.second += y.second;
  return x;
}

template <class f>
inline auto& operator-=(R<f>& x, f y) {
  x.first -= y;
  return x;
}

template <class f>
inline auto operator*(R<f> x, R<f> y) {
  return R<f>{x.first * y.first,
              x.second * (y.first > 0 ? y.first : -y.first) +
                  y.second * (x.first > 0 ? x.first : -x.first)};
}

template <class f>
inline auto operator*(R<f> x, f y) {
  return R<f>{x.first * y, x.second * (y > 0 ? y : -y)};
}

template <class f>
inline auto operator*(f x, R<f> y) {
  return R<f>{x * y.first, y.second * (x > 0 ? x : -x)};
}

template <class f>
inline auto& operator*=(R<f>& x, R<f> y) {
  x.second = x.second * (y.first > 0 ? y.first : -y.first) +
             y.second * (x.first > 0 ? x.first : -x.first);
  x.first *= y.first;
  return x;
};

template <class f>
inline auto& operator*=(R<f>& x, f y) {
  x.second = x.second * (y > 0 ? y : -y);
  x.first *= y;
  return x;
};

template <class f>
inline auto operator/(R<f> x, R<f> y) {
  return R<f>{x.first / y.first, x.second / (y.first > 0 ? y.first : -y.first) +
                                     (x.first > 0 ? x.first : -x.first) *
                                         y.second / (y.first * y.first)};
}

template <class f>
inline auto operator/(R<f> x, f y) {
  return R<f>{x.first / y, x.second / (y > 0 ? y : -y)};
}

template <class f>
inline auto operator/(f x, R<f> y) {
  return R<f>{x / y.first, (x > 0 ? x : -x) * y.second / (y.first * y.first)};
}

template <class f>
inline auto& operator/=(R<f>& x, R<f> y) {
  x.second =
      x.second / (y.first > 0 ? y.first : -y.first) +
      (x.first > 0 ? x.first : -x.first) * y.second / (y.first * y.first);
  x.first /= y.first;
  return x;
}  // namespace ffd::R

template <class f>
inline auto& operator/=(R<f>& x, f y) {
  x.second = x.second / (y > 0 ? y : -y);
  x.first /= y;
  return x;
}  // namespace ffd::R

template <class f>
inline auto operator+(R<f> x) {
  return x;
}

template <class f>
inline auto operator-(R<f> x) {
  return R<f>{-x.first, x.second};
}

template <class f>
inline auto abs(R<f> x) {
  return (x.first > 0 ? x.first : -x.first);
}

}  // namespace ffd::round_off
