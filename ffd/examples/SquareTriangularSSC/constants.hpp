template <int>
struct int_type {
  int j;
  bool constexpr operator==(int_type x) const { return x.j == j; }
  bool constexpr operator!=(int_type x) const { return x.j != j; }
};
template <int>
struct char_type {
  char j;
  operator char() const { return j; }
  bool constexpr operator==(char_type x) const { return x.j == j; }
  bool constexpr operator!=(char_type x) const { return x.j != j; }
};

char_type<0> constexpr square{'s'};
char_type<0> constexpr triangle{'t'};
auto constexpr triangular = triangle;
auto constexpr rectangle = square;
auto constexpr rectangular = square;
