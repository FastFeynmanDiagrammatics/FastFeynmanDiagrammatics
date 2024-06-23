namespace ffd::canonical_sigma::unit_test {

auto lexicographical_permutation() {
  unsigned int v;  // current permutation of bits
  unsigned int w;  // next permutation of bits
  int n = 6;
  int k = 1;
  for (ulong k = 1; k <= n; ++k) {
    int counter = 0;
    for (ulong S = (1ul << k) - 1; S < (1ul << n);
         S = ffd::set_theory::next_permutation(S)) {
      counter++;
      std::cerr << k << " " << S << " " << __builtin_popcount(S) << "\n";
    }  // for S in range((1ul<<k)-1, (1ul<<n))
    std::cerr << "counter_" << k << " = " << counter << "\n";
  }
}

}  // namespace ffd::canonical_sigma::unit_test
