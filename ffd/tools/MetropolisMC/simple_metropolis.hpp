namespace ffd::metropolis_mc {

  template <class Field_f = Real>
  auto step(
	    std::array<Field_f, 2> f
	    ) {
    using std::abs;
    // perform metropolis update with some probability
    const auto ratio = abs(f[1]) / abs(f[0]);
    auto c = (ratio > random_distributions::Proba());
    return c;
  }


} // namespace ffd
