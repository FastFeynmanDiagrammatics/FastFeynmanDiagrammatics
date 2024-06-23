namespace ffd::alpha_f {

template <typename PM_t, PM_shift_t>
auto AddDiagonalPM(PM_t const& PM, PM_shift_t& PM_shift, std::size_t n) {
  auto const size_PM = PM.size();

  for (uint S = 0; S < size_PM; ++S) {
    PM_shift[S][0] = PM[S + shift];
  }
  for (uint j = 1; j < card; j *= 2) {
    for (uint S = card - 1; S >= j; --S) {
      if (j & S) {
        const auto alt_S = S - j;
        for (uint k = 1; k <= std::min(__builtin_popcount(S), half_n); ++k) {
          PM_shift[S][k] += PM_shift[alt_S][k - 1];
        }
      }
    }
  }
  return P_array;
}

}  // namespace ffd::alpha_f
