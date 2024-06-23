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

char constexpr square = 's';
char constexpr triangle = 't';
auto constexpr triangular = triangle;

char constexpr BZ = 'z';
char constexpr QQ = 'q';
char constexpr QP = 'p';
char constexpr Q0 = 'o';
