namespace ffd::random_distributions {

// usage: to get a number out of these probability distributions,
// one uses Distribution(RNGen)
std::random_device RandomDevice1;
std::mt19937_64 RNGen(RandomDevice1());
std::uniform_int_distribution<int> ZeroOr1Distribution(0, 1);
std::uniform_real_distribution<Real> Proba0_1(0, 1);
std::normal_distribution<Real> GaussianDistribution(0, 1);

inline Real Proba() {
  return Proba0_1(RNGen);
}

inline bool FlipCoin() {
  return Proba() > .5;
}

inline int RandomSign() {
  return 2 * FlipCoin() - 1;
}

// one integer from j1 to j2 (both included)
inline int RandomIntFrom1To2(int j1, int j2) {
  assert((j1 <= j2));
  std::uniform_int_distribution<int> dist_temp(j1, j2);
  return dist_temp(RNGen);
}

// one integer from 0 to j (both included)
inline int RandomIntFromZeroTo(int j) {
  return RandomIntFrom1To2(0, j);
}

template <typename IntType1, typename IntType2>
inline IntType1 RandomInRange(IntType1 j1, IntType2 j2) {
  assert((j1 < IntType1(j2)));
  std::uniform_int_distribution<IntType1> dist_temp(j1, j2 - 1);
  return dist_temp(RNGen);
}

template <typename IntType>
inline IntType RandomInRange(IntType j) {
  return RandomInRange<IntType, IntType>(0, j);
}

std::string RandomString(std::size_t length) {
  auto randchar = []() -> char {
    const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const std::size_t max_index = (sizeof(charset));
    return charset[RandomInRange(max_index - 1)];
  };
  std::string str(length, 0);
  std::generate_n(str.begin(), length, randchar);
  return str;
}

}  // namespace ffd::random_distributions

namespace ffd::testing_space {

using ffd::random_distributions::Proba;

using ffd::random_distributions::RNGen;

}  // namespace ffd::testing_space
