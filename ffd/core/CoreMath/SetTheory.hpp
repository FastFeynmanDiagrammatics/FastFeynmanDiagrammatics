namespace ffd::set_theory {

constexpr int CardinalitySet(BinaryInt S) {
  BinaryInt ret = 0;
  for (; S != 0; ++ret, S &= S - 1) {
  }
  return ret;
}

auto VectorOfBinaryDigitsOf(BinaryInt S) {
  std::vector<BinaryInt> ret(CardinalitySet(S), 0u);
  for (BinaryInt j = 0, counter = 0; S > 0; S >>= 1, ++j) {
    ret[counter] = j;
    counter += S & 1u;
  }
  return ret;
}

template <class vector_t>
void BinaryDigitsOf_r(ulong S, vector_t& vec) {
  for (ulong j = 0, counter = 0; S > 0; S >>= 1, ++j) {
    vec[counter] = j;
    counter += S & 1u;
  }
}

uint Order(uint v) {
  uint r = (v > 0);
  while (v >>= 1) {
    r++;
  }
  return r;
}

uint BinarySpread(uint set, uint step) {
  uint val = 0;
  for (uint v = 0; v <= Order(set); v++) {
    if (set & (1 << v)) {
      for (uint st = 0; st < step; ++st) {
        val += (1 << (step * (v + 1) - 1 - st));
      }
    }
  }
  return val;
}

inline unsigned next_permutation(uint const v) {
  uint const t = v | (v - 1);
  return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
}

}  // namespace ffd::set_theory
