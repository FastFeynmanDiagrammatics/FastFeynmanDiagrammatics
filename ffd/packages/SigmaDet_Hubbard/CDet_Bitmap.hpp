namespace ffd::sigmadet_hubbard {

template <std::size_t n>
auto CDet_Bitmap() {
  array2d<std::uint8_t, n + 1, (1ul << n)> ret;
  ret.fill(0);
  for (ulong S = 0; S < (1ul << n); ++S) {
    std::vector<BinaryInt> const digits_S =
        ffd::set_theory::VectorOfBinaryDigitsOf(S);
    ret(0, S) = __builtin_popcount(S);
    for (ulong j = 0; j < size(digits_S); ++j) {
      ret(j + 1, S) = digits_S[j];
    }  // for j in range(0, n)
  }    // for S in range(0, (1yl<<n))
  return ret;
}

}  // namespace ffd::sigmadet_hubbard
