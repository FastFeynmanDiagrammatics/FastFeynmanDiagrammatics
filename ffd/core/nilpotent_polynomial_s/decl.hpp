namespace ffd::nilpotent_polynomial_s {

enum mode { unsafe, safe };

template <class Field = Real, auto n = 0, mode m = safe>
struct P {
  using field = Field;
  Real precision = std::numeric_limits<field>::epsilon() * std::pow(10., 5.);

  std::array<field, (1ul << n)> coef;

  inline std::size_t size() const { return coef.size(); }
  inline auto constexpr order() const { return n; }

  P() {}

  P(field x) {
    coef[0] = x;
    coef.fill(field(0));
  }

  inline field& operator[](BinaryInt S) { return coef[S]; }
  inline field const& operator[](BinaryInt S) const { return coef[S]; }
};

}  // namespace ffd::nilpotent_polynomial_s
